// amrlpt_wrapper.cpp
// C++ bridge between Fortran amrlpt_class and AMReX NeighborParticleContainer
// Particle struct: 16 extra reals (d, vel[3], angVel[3], Acol[3], Tcol[3], dt, mass, eint)
//                  1 extra int   (flag)
// All functions callable from Fortran via bind(C)
//
// Uses NeighborParticleContainer (backward-compatible superset of
// AmrParticleContainer):
//   - fillNeighbors / clearNeighbors: ghost particles for collision detection
//   - buildNeighborList: explicit pair lists for DEM / peridynamics
//
// Ghost particle layout (after fillNeighbors):
//   ptile.GetArrayOfStructs()[0 .. numRealParticles()-1]         : valid
//   particles ptile.GetArrayOfStructs()[numRealParticles() .. total-1]     :
//   ghost particles
//
// Neighbor list (after buildNeighborList):
//   CSR format: m_nbor_offsets[np+1], m_nbor_list[ntotal], both unsigned int
//   For valid particle i, neighbors are m_nbor_list[offsets[i] ..
//   offsets[i+1]-1] Index j < numRealParticles() => valid particle; j >=
//   numRealParticles() => ghost

#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_NeighborList.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <cstddef>
#include <type_traits>

using namespace amrex;

// -----------------------------------------------------------------------
// Particle container type: 16 extra reals, 1 extra int
// Full in-memory layout per particle (double precision):
//   pos[3]     (24 B)  -- AMReX base
//   rdata[16]  (128 B) -- d, vel[3], angVel[3], Acol[3], Tcol[3], dt, mass, eint
//   idcpu      (8 B)   -- AMReX packed id+cpu
//   idata[1]   (4 B)   -- flag
// -----------------------------------------------------------------------
#define AMRLPT_NREAL 16
#define AMRLPT_NINT 1

namespace {

// -----------------------------------------------------------------------
//   AMRLPTPC: Thin subclass of NeighborParticleContainer
// Purposes:
//   1. Forward-declare AmrCore-based constructor (ParGDB pass-through)
//   2. Expose setNeighborCells to allow per-call ghost radius changes
//   3. Expose m_neighbor_list per tile (for CSR neighbor list access)
// -----------------------------------------------------------------------
class AMRLPTPC : public NeighborParticleContainer<AMRLPT_NREAL, AMRLPT_NINT> {
public:
  using Base = NeighborParticleContainer<AMRLPT_NREAL, AMRLPT_NINT>;
  using ParticleType = Base::ParticleType;

  explicit AMRLPTPC(AmrCore *amrcore)
      : Base(amrcore->GetParGDB(),
             1) // ncells=1 default; overridden via setNeighborCells
  {}

  // Allow changing the ghost radius before each fillNeighbors() call
  void setNeighborCells(int n) { m_num_neighbor_cells = n; }

  // CSR access to the neighbor list for a tile (built by buildNeighborList)
  NeighborList<ParticleType> &getNeighborList(int lev, int grid, int tile) {
    return m_neighbor_list[lev][std::make_pair(grid, tile)];
  }
};

using PC = AMRLPTPC;
using PT = PC::ParticleType;
static_assert(AMREX_SPACEDIM == 3 && sizeof(ParticleReal) == 8,
              "Fortran particle bridge requires 3D, double-precision particles");
static_assert(std::is_standard_layout<PT>::value && sizeof(PT) == 168,
              "Particle<16,1> no longer matches the Fortran bind(C) layout");
static_assert(offsetof(PT, m_rdata) == 24 && offsetof(PT, m_idcpu) == 152 &&
              offsetof(PT, m_idata) == 160, "Particle member offsets changed");
constexpr const char* layout_file = "NGA2ParticleLayout";
constexpr const char* layout_tag = "NGA2_AMRLPT_MASS_EINT";
constexpr int layout_version = 1;

} // namespace

extern "C" {

std::size_t amrlpt_particle_size() { return sizeof(PT); }

// -----------------------------------------------------------------------
// Lifecycle
// -----------------------------------------------------------------------

void amrlpt_new_pc(PC *&pc, void *amrcore_raw) {
  pc = new PC(static_cast<AmrCore *>(amrcore_raw));
}

void amrlpt_delete_pc(PC *pc) { delete pc; }

// -----------------------------------------------------------------------
// Redistribution (AMR-aware particle sorting to finest covering level)
// -----------------------------------------------------------------------

void amrlpt_redistribute(PC *pc, int lev_min, int lev_max, int ng) {
  pc->Redistribute(lev_min, lev_max, ng);
}

// Enforce an NGA2 particle-level ceiling even if particles already reside above
// it.
void amrlpt_enforce_max_level(PC *pc, int lev_max, long long &number_demoted) {
  if (pc == nullptr) {
    amrex::Abort("amrlpt_enforce_max_level: null particle container");
  }

  const int finest = pc->finestLevel();
  if (lev_max < 0 || lev_max > finest) {
    amrex::Abort(
        "amrlpt_enforce_max_level: level cap is outside the AMR hierarchy");
  }

  // Neighbor particles are communication copies appended after the real
  // particles.  Remove them before copying or changing any particle tile.
  pc->clearNeighbors();

  const Long count_before = pc->TotalNumberOfParticles(true, false);
  Long local_demoted = 0;
  auto &levels = pc->GetParticles();
  const int stored_finest = static_cast<int>(levels.size()) - 1;

  // Count first so the staging tile is allocated only once.
  for (int lev = lev_max + 1; lev <= stored_finest; ++lev) {
    for (auto const &kv : levels[lev]) {
      local_demoted += static_cast<Long>(kv.second.numRealParticles());
    }
  }
  if (local_demoted > static_cast<Long>(std::numeric_limits<int>::max())) {
    amrex::Abort("amrlpt_enforce_max_level: too many local particles for one "
                 "staging tile");
  }

  PC::ParticleTileType staged;
  staged.resize(static_cast<int>(local_demoted));
  auto &staged_aos = staged.GetArrayOfStructs();
  Long offset = 0;

  // Copy only real particles.  PT assignment preserves the AMReX ID and all
  // structure data; unlike normal particle injection, no new ID is assigned.
  for (int lev = lev_max + 1; lev <= stored_finest; ++lev) {
    for (auto const &kv : levels[lev]) {
      auto const &src_tile = kv.second;
      auto const &src_aos = src_tile.GetArrayOfStructs();
      const int np = src_tile.numRealParticles();
      for (int i = 0; i < np; ++i) {
        staged_aos[static_cast<int>(offset++)] = src_aos[i];
      }
    }
  }

  for (int lev = lev_max + 1; lev <= stored_finest; ++lev) {
    pc->RemoveParticlesAtLevel(lev);
  }

  // AddParticlesAtLevel is collective and is safe with an empty staging tile.
  // Its internal level-zero redistribution gives every staged particle a
  // valid tile before the capped hierarchy redistribution below.
  pc->AddParticlesAtLevel(staged, 0);
  pc->Redistribute(0, lev_max, 0);

  const Long count_after = pc->TotalNumberOfParticles(true, false);
  if (count_after != count_before) {
    amrex::Abort("amrlpt_enforce_max_level: particle count changed while "
                 "enforcing level cap");
  }

  Long remaining_above = 0;
  for (int lev = lev_max + 1; lev <= stored_finest; ++lev) {
    for (auto const &kv : pc->GetParticles(lev)) {
      remaining_above += static_cast<Long>(kv.second.numRealParticles());
    }
  }
  ParallelAllReduce::Sum(local_demoted, ParallelContext::CommunicatorSub());
  ParallelAllReduce::Sum(remaining_above, ParallelContext::CommunicatorSub());
  if (remaining_above != 0) {
    amrex::Abort(
        "amrlpt_enforce_max_level: particles remain above the requested cap");
  }
  if (!pc->OK(0, lev_max, 0)) {
    amrex::Abort("amrlpt_enforce_max_level: capped particle container failed "
                 "AMReX validity check");
  }

  number_demoted = static_cast<long long>(local_demoted);
}

// Global number of real particles stored on each requested AMR level.
void amrlpt_level_counts(PC *pc, long long *counts, int nlevels) {
  if (pc == nullptr || counts == nullptr || nlevels < 0) {
    amrex::Abort("amrlpt_level_counts: invalid argument");
  }

  Vector<Long> work(static_cast<std::size_t>(nlevels), Long(0));
  const int stored_levels = static_cast<int>(pc->GetParticles().size());
  const int nscan = std::min(nlevels, stored_levels);
  for (int lev = 0; lev < nscan; ++lev) {
    for (auto const &kv : pc->GetParticles(lev)) {
      work[lev] += static_cast<Long>(kv.second.numRealParticles());
    }
  }
  if (nlevels > 0) {
    ParallelAllReduce::Sum(work.data(), nlevels,
                           ParallelContext::CommunicatorSub());
  }
  for (int lev = 0; lev < nlevels; ++lev) {
    counts[lev] = static_cast<long long>(work[lev]);
  }
}

// -----------------------------------------------------------------------
// Ghost (neighbor) particles for collision detection and short-range
// interactions
//
// fillNeighbors communicates particles within ngrow cells of each tile boundary
// into the tile's particle array (appended after numRealParticles).
// fillNeighborsRadius only communicates particles within physical distance
// search_radius of tile boundaries (much less data when search_radius << dx).
// clearNeighbors removes them.
// -----------------------------------------------------------------------

void amrlpt_fill_neighbors(PC *pc, int ngrow) {
  pc->setNeighborCells(ngrow);
  pc->fillNeighbors();
}

void amrlpt_fill_neighbors_radius(PC *pc, double search_radius) {
  pc->fillNeighbors(amrex::Real(search_radius));
}

void amrlpt_clear_neighbors(PC *pc) { pc->clearNeighbors(); }

// -----------------------------------------------------------------------
// Ghost particle access per tile.
// Call after amrlpt_fill_neighbors.
// Ghost particles are appended at ptile[numRealParticles() .. numTotal-1].
// -----------------------------------------------------------------------

void amrlpt_get_neighbor_particles_mfi(PC *pc, int lev, MFIter *mfi, PT *&dp,
                                       long long &np) {
  if (lev >= static_cast<int>(pc->GetParticles().size())) {
    np = 0;
    dp = nullptr;
    return;
  }
  const int grid = mfi->index();
  const int tile = mfi->LocalTileIndex();
  auto &plev = pc->GetParticles(lev);
  auto it = plev.find(std::make_pair(grid, tile));
  if (it != plev.end()) {
    auto &ptile = it->second;
    np = static_cast<long long>(ptile.numNeighborParticles());
    dp = (np > 0) ? ptile.GetArrayOfStructs().data() + ptile.numRealParticles()
                  : nullptr;
  } else {
    np = 0;
    dp = nullptr;
  }
}

// -----------------------------------------------------------------------
// Neighbor list (explicit pair list) for DEM / peridynamics.
// buildNeighborList: must be called AFTER fillNeighbors.
//   check_pair criterion: |r_i - r_j| < rcrit
// Neighbor list is per-tile in CSR format:
//   offsets[0..np]: prefix sums (unsigned int)
//   list[0..ntot-1]: neighbor indices into full particle array (valid+ghost)
//   np: number of valid particles
//   ntot: total neighbor entries
// -----------------------------------------------------------------------

void amrlpt_build_neighbor_list(PC *pc, double rcrit) {
  const double rcrit2 = rcrit * rcrit;
  auto check_pair = [rcrit2](const PT &p1, const PT &p2) -> bool {
    double d2 = 0.0;
    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
      d2 += (p1.pos(dim) - p2.pos(dim)) * (p1.pos(dim) - p2.pos(dim));
    return d2 < rcrit2;
  };
  pc->buildNeighborList(check_pair, amrex::Real(rcrit), false);
}

// Returns CSR offsets and list for the neighbor list per tile.
// offsets: pointer to np+1 unsigned ints (0-indexed relative offsets into list)
// list:    pointer to ntot unsigned ints (indices into valid+ghost particle
// array) np:      number of valid particles (list driven from particle i in
// [0,np-1]) ntot:    total number of neighbor entries across all particles
void amrlpt_get_neighbor_list_mfi(PC *pc, int lev, MFIter *mfi,
                                  const unsigned int *&offsets,
                                  const unsigned int *&list, long long &np,
                                  long long &ntot) {
  const int grid = mfi->index();
  const int tile = mfi->LocalTileIndex();
  auto &nl = pc->getNeighborList(lev, grid, tile);
  np = static_cast<long long>(nl.numParticles());
  ntot = static_cast<long long>(nl.GetList().size());
  offsets = (np > 0) ? nl.GetOffsets().dataPtr() : nullptr;
  list = (ntot > 0) ? nl.GetList().dataPtr() : nullptr;
}

// -----------------------------------------------------------------------
// Particle access via MFIter (grid+tile indices come from the MFIter)
// Returns pointer to VALID particle array + count for this tile.
// Ghost particles (if present after fillNeighbors) are excluded here;
// use amrlpt_get_neighbor_particles_mfi for them.
// -----------------------------------------------------------------------

void amrlpt_get_particles_mfi(PC *pc, int lev, MFIter *mfi, PT *&dp,
                              long long &np) {
  // Guard: m_particles is empty if Redistribute has never been called
  if (lev >= static_cast<int>(pc->GetParticles().size())) {
    np = 0;
    dp = nullptr;
    return;
  }
  const int grid = mfi->index();
  const int tile = mfi->LocalTileIndex();
  auto &plev = pc->GetParticles(lev);
  auto it = plev.find(std::make_pair(grid, tile));
  if (it != plev.end()) {
    auto &ptile = it->second;
    np = static_cast<long long>(ptile.numRealParticles()); // excludes ghosts
    dp = (np > 0) ? ptile.GetArrayOfStructs().data() : nullptr;
  } else {
    np = 0;
    dp = nullptr;
  }
}

// Returns the full combined (valid+ghost) particle array for a tile,
// along with np_total (valid+ghost count) and np_valid (valid count only).
// np_valid is the count to use for the outer initiating loop;
// np_total covers the full index space referenced by the neighbor list.
void amrlpt_get_all_particles_mfi(PC *pc, int lev, MFIter *mfi, PT *&dp,
                                  long long &np_total, long long &np_valid) {
  // Guard: m_particles is empty if Redistribute has never been called
  if (lev >= static_cast<int>(pc->GetParticles().size())) {
    np_total = 0;
    np_valid = 0;
    dp = nullptr;
    return;
  }
  const int grid = mfi->index();
  const int tile = mfi->LocalTileIndex();
  auto &plev = pc->GetParticles(lev);
  auto it = plev.find(std::make_pair(grid, tile));
  if (it != plev.end()) {
    auto &ptile = it->second;
    np_valid = static_cast<long long>(ptile.numRealParticles());
    np_total =
        static_cast<long long>(ptile.numTotalParticles()); // real + neighbor
    dp = (np_total > 0) ? ptile.GetArrayOfStructs().data() : nullptr;
  } else {
    np_total = 0;
    np_valid = 0;
    dp = nullptr;
  }
}

void amrlpt_num_particles_mfi(PC *pc, int lev, MFIter *mfi, long long &np) {
  if (lev >= static_cast<int>(pc->GetParticles().size())) {
    np = 0;
    return;
  }
  const int grid = mfi->index();
  const int tile = mfi->LocalTileIndex();
  auto &plev = pc->GetParticles(lev);
  auto it = plev.find(std::make_pair(grid, tile));
  np = (it != plev.end())
           ? static_cast<long long>(it->second.numRealParticles())
           : 0;
}

// -----------------------------------------------------------------------
// Particle grid access (ParGDB: particle-specific BA/DM)
// In single-grid mode these fall back to the fluid grid.
// In dual-grid mode these return the particle-specific grid.
// -----------------------------------------------------------------------

void amrlpt_get_particle_boxarray(PC *pc, int lev, void **ba_ptr) {
  *ba_ptr = const_cast<amrex::BoxArray *>(&(pc->ParticleBoxArray(lev)));
}

void amrlpt_get_particle_distromap(PC *pc, int lev, void **dm_ptr) {
  *dm_ptr = const_cast<amrex::DistributionMapping *>(
      &(pc->ParticleDistributionMap(lev)));
}

void amrlpt_set_particle_boxarray(PC *pc, int lev, void *ba_ptr) {
  auto *ba = static_cast<amrex::BoxArray *>(ba_ptr);
  pc->SetParticleBoxArray(lev, *ba);
}

void amrlpt_set_particle_distromap(PC *pc, int lev, void *dm_ptr) {
  auto *dm = static_cast<amrex::DistributionMapping *>(dm_ptr);
  pc->SetParticleDistributionMap(lev, *dm);
}

// -----------------------------------------------------------------------
// Add a single particle to a tile (by explicit grid+tile index)
// -----------------------------------------------------------------------

void amrlpt_add_particle_i(PC *pc, int lev, int grid, int tile, PT *p) {
  if (lev >= static_cast<int>(pc->GetParticles().size()))
    return;
  auto &plev = pc->GetParticles(lev);
  plev[std::make_pair(grid, tile)].push_back(*p);
}

// -----------------------------------------------------------------------
// Global unique particle ID counter and current CPU rank
// -----------------------------------------------------------------------

void amrlpt_get_next_id(long long &id) { id = PT::NextID(); }

void amrlpt_set_next_id(long long id) { PT::NextID(id); }

void amrlpt_get_cpu(int &cpu) { cpu = ParallelDescriptor::MyProc(); }

// -----------------------------------------------------------------------
// ID/CPU accessors for the packed idcpu field (private in Fortran struct)
// -----------------------------------------------------------------------

void amrlpt_get_particle_id(long long &id, const PT *p) { id = p->id(); }

void amrlpt_set_particle_id(long long id, PT *p) { p->id() = id; }

void amrlpt_get_particle_cpu(int &cpu, const PT *p) { cpu = p->cpu(); }

void amrlpt_set_particle_cpu(int cpu, PT *p) { p->cpu() = cpu; }

void amrlpt_particle_is_valid(int &valid, const PT *p) {
  valid = p->id().is_valid() ? 1 : 0;
}

// -----------------------------------------------------------------------
// Total particle count (global, across all ranks and levels)
// -----------------------------------------------------------------------

void amrlpt_total_np(PC *pc, long long &np) {
  np = pc->TotalNumberOfParticles();
}

// -----------------------------------------------------------------------
// Checkpoint: write particles to fullpath (e.g. "restart/part_1.00E+00").
// Splits at last '/' to satisfy AMReX Checkpoint(dir, name) API.
// Creates parent directory hierarchy via PreBuildDirectorHierarchy.
// -----------------------------------------------------------------------

void amrlpt_write(PC *pc, const char *fullpath, int is_chk) {
  std::string path(fullpath);
  while (!path.empty() && path.back() == '/')
    path.pop_back();
  auto pos = path.rfind('/');
  std::string parent =
      (pos != std::string::npos) ? path.substr(0, pos) : std::string(".");
  std::string leaf = (pos != std::string::npos) ? path.substr(pos + 1) : path;
  amrex::PreBuildDirectorHierarchy(parent, "", 1, true);
  pc->Checkpoint(parent, leaf, is_chk != 0);
  if (ParallelDescriptor::IOProcessor()) {
    std::ofstream schema(path + "/" + layout_file);
    schema << layout_tag << ' ' << layout_version << ' ' << AMRLPT_NREAL
           << ' ' << AMRLPT_NINT << '\n';
    schema.close();
    if (!schema) amrex::Abort("amrlpt_write: failed to write particle layout metadata");
  }
  ParallelDescriptor::Barrier();
}

// -----------------------------------------------------------------------
// Visualization plotfile with selective component output.
// write_real[AMRLPT_NREAL]: bitmask (1=write, 0=skip) for each extra real.
// write_int[AMRLPT_NINT]:   bitmask for each extra int.
// Component names are hardcoded here to match the Fortran part struct layout:
//   rdata[0..15]: d, vx, vy, vz, wx, wy, wz, ax, ay, az, tx, ty, tz, dt, mass, eint
//   idata[0]:     flag
// Position (pos[3]) is always written by AMReX (baked into particle format).
// -----------------------------------------------------------------------

void amrlpt_write_plotfile(PC *pc, const char *basedir, const char *pname,
                           const int *write_real, const int *write_int,
                           double time) {
  static const Vector<std::string> rnames = {"d",  "vx", "vy", "vz", "wx",
                                             "wy", "wz", "ax", "ay", "az",
                                             "tx", "ty", "tz", "dt", "mass", "eint"};
  static const Vector<std::string> inames = {"flag"};

  Vector<int> wr(write_real, write_real + AMRLPT_NREAL);
  Vector<int> wi(write_int, write_int + AMRLPT_NINT);

  pc->WritePlotFile(std::string(basedir), std::string(pname), wr, wi, rnames,
                    inames);

  // Store simulation time inside the plotfile directory (IOProcessor only)
  if (ParallelDescriptor::IOProcessor()) {
    std::ofstream tf(std::string(basedir) + "/time");
    tf << std::setprecision(17) << time << '\n';
  }
}

// -----------------------------------------------------------------------
// Read back the simulation time stored by amrlpt_write_plotfile.
// Returns -1.0 if the file does not exist or cannot be read.
// -----------------------------------------------------------------------

double amrlpt_read_plotfile_time(const char *basedir) {
  std::ifstream tf(std::string(basedir) + "/time");
  double t = -1.0;
  if (tf.good())
    tf >> t;
  return t;
}

// -----------------------------------------------------------------------
// Restart: read back a particle checkpoint written by amrlpt_write.
// -----------------------------------------------------------------------

void amrlpt_read(PC *pc, const char *fullpath) {
  // Mirror amrlpt_write: split a single fullpath into (parent, leaf) so the
  // ABI matches the Fortran caller, which passes one path argument.
  std::string path(fullpath);
  while (!path.empty() && path.back() == '/')
    path.pop_back();
  auto pos = path.rfind('/');
  std::string parent =
      (pos != std::string::npos) ? path.substr(0, pos) : std::string(".");
  std::string leaf = (pos != std::string::npos) ? path.substr(pos + 1) : path;
  // Reject pre-mass/eint checkpoints before AMReX interprets their particle payload.
  std::ifstream schema(path + "/" + layout_file);
  std::string tag;
  int version = 0, nr = 0, ni = 0;
  if (!(schema >> tag >> version >> nr >> ni) || tag != layout_tag ||
      version != layout_version || nr != AMRLPT_NREAL || ni != AMRLPT_NINT) {
    amrex::Abort("amrlpt_read: incompatible particle checkpoint layout; expected "
                 "MASS_EINT v1 (16 reals, 1 integer). Legacy 14-real checkpoints "
                 "require explicit migration; start a fresh run.");
  }
  // Restart replaces the current state, including when called on a populated container.
  pc->clearParticles();
  pc->Restart(parent, leaf);
}

// -----------------------------------------------------------------------
// Append an array of new particles to the container at level 0.
// Called collectively by ALL ranks.  Any rank may pass n>0 with valid data
// Ranks with nothing to add pass n=0 and raw=nullptr.
// -----------------------------------------------------------------------

void amrlpt_append_particles(PC *pc, const void *raw, long long n) {
  PC::ParticleTileType ptile;

  if (n > 0 && raw != nullptr) {
    const PT *src = reinterpret_cast<const PT *>(raw);
    ptile.resize(static_cast<int>(n));
    auto &aos = ptile.GetArrayOfStructs();
    for (long long i = 0; i < n; ++i) {
      PT p = src[i];         // copy pos + rdata + idata
      p.id() = PT::NextID(); // assign valid AMReX identity
      p.cpu() = ParallelDescriptor::MyProc();
      aos[static_cast<int>(i)] = p;
    }
  }

  pc->AddParticlesAtLevel(ptile, 0); // collective: level 0 redistribute inside
}

} // extern "C"
