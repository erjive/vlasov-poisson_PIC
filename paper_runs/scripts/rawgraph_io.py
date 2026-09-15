"""
Reader for VP_PIC's raw binary output (output_format="raw",
"<directory>/vlasov_output.raw", written by src/raw_io.f90).

No self-description (unlike HDF5) -- this module hardcodes the exact
byte layout documented in raw_io.f90's header comment, and must be kept
in sync with it by hand if that layout ever changes.

Layout (all int64/float64, native machine endianness -- little-endian
on x86_64):

  HEADER (once):     Nr:int64, autointeraction:int64(0/1), r:float64[Nr]
  RECORD (per step):  l:int64, time,kinetic,potential,total_energy:float64,
                       Npart:int64, rho,avg_rho,curr:float64[Nr],
                       [force,potential:float64[Nr] if autointeraction],
                       r_part,p_part,f:float64[Npart]

Records are NOT fixed-size (Npart can shrink via reduceparticles), so
the file is indexed by one sequential pass reading just the int64
Npart field of each record and seeking over the rest -- cheap even for
thousands of records, since it never touches the bulk float64 payload
until a specific record is actually requested.
"""
import numpy as np


class RawRun:
    """Indexes a vlasov_output.raw file and reads snapshots on demand,
    lazily (like h5py) -- nothing but the header and a byte-offset index
    is held in memory until you ask for a specific record's arrays."""

    def __init__(self, path):
        self.path = path
        self._fh = open(path, "rb")
        self.Nr = int(np.fromfile(self._fh, dtype=np.int64, count=1)[0])
        self.autointeraction = bool(np.fromfile(self._fh, dtype=np.int64, count=1)[0])
        self.r = np.fromfile(self._fh, dtype=np.float64, count=self.Nr)
        self._grid_stat_names = ("rho", "avg_rho", "curr") + (
            ("force", "potential") if self.autointeraction else ()
        )
        self._index = []  # list of dicts: offset, l, time, energies, Npart
        self._build_index()

    def _record_nbytes(self, npart):
        header = 8 * (1 + 4 + 1)  # l, time, kinetic, potential, total_energy, Npart
        grid = 8 * self.Nr * len(self._grid_stat_names)
        parts = 8 * npart * 3  # r_part, p_part, f
        return header + grid + parts

    def _build_index(self):
        fh = self._fh
        while True:
            start = fh.tell()
            head = np.fromfile(fh, dtype=np.int64, count=1)
            if head.size == 0:
                break
            l = int(head[0])
            time, kinetic, potential, total_energy = np.fromfile(fh, dtype=np.float64, count=4)
            npart = int(np.fromfile(fh, dtype=np.int64, count=1)[0])
            self._index.append(dict(offset=start, l=l, time=float(time),
                                     kinetic=float(kinetic), potential=float(potential),
                                     total_energy=float(total_energy), npart=npart))
            payload = 8 * self.Nr * len(self._grid_stat_names) + 8 * npart * 3
            fh.seek(payload, 1)

    @property
    def n_steps(self):
        return len(self._index)

    @property
    def times(self):
        return [e["time"] for e in self._index]

    def read(self, i):
        """Read snapshot i (0-based, in save order) as a dict of arrays
        plus the scalar attributes, mirroring the h5py group interface
        used elsewhere (g['rho'][:], g.attrs['time'])."""
        e = self._index[i]
        fh = self._fh
        fh.seek(e["offset"] + 8 * 6)  # skip l,time,kinetic,potential,total_energy,Npart
        out = {"attrs": {k: e[k] for k in ("time", "kinetic", "potential", "total_energy", "l")}}
        for name in self._grid_stat_names:
            out[name] = np.fromfile(fh, dtype=np.float64, count=self.Nr)
        npart = e["npart"]
        out["r_part"] = np.fromfile(fh, dtype=np.float64, count=npart)
        out["p_part"] = np.fromfile(fh, dtype=np.float64, count=npart)
        out["f"] = np.fromfile(fh, dtype=np.float64, count=npart)
        return out

    def close(self):
        self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()
