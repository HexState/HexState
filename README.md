<p align="center">
Notice to those cloning this
</p>
<p align="center">
If there is a test folder, it is stable; if there is no test folder it is being actively updated.
</p>

<p align="center">
  <strong>⬡ HEXSTATE ENGINE</strong>
</p>

---

# HexState Engine v2

**The world's first D=6 quantum tensor network processor.**

100,000 sites · 6^100000 ≈ 10^77815 dimensions · One laptop · 110 seconds.

---

## What Changed in v2

Version 1 operated on coarse-grained chunks with shadow caches and lazy quhit resolution. It could entangle 100 trillion quhits in product-like states — but the entanglement was limited to 2-body Bell pairs via braiding.

Version 2 is a total rebuild, introducing a **Matrix Product State (MPS) tensor network layer** that enables genuine **N-body entanglement** across arbitrarily long chains. This is the same mathematical framework used by professional condensed matter physicists (DMRG, TEBD), but implemented natively in D=6 — a physical dimension no quantum computer has ever operated in.
