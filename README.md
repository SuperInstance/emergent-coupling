# emergent-coupling

**Spectral gap emergence — when two systems couple, does their joint structure exceed the sum of parts?**

Investigates emergent behavior when independent spectral systems (each with eigenvalues and dominant modes) are coupled. Computes whether the coupled system's spectral gap, conservation ratio, and synergy exceed what either system achieves alone.

## What This Gives You

- **Coupling mode detection** — classifies coupled systems as Reinforcing, Complementary, or Competitive
- **Emergence detection** — flags when coupled spectral gap exceeds both individual gaps
- **Synergy score** — quantitative measure of spectral coupling benefit
- **Alignment metric** — cosine similarity between dominant modes
- **Conservation ratio tracking** — λ₂/λₙ for individual and coupled systems
- **Zero dependencies** — pure Rust

## Quick Start

```rust
use emergent_coupling::{SpectralSystem, couple};

let sys_a = SpectralSystem::new("oscillator", 
    vec![0.0, 0.3, 1.5, 3.0],  // eigenvalues
    vec![0.5, 0.5, 0.3, -0.7],  // dominant eigenvector
);

let sys_b = SpectralSystem::new("rotor",
    vec![0.0, 0.1, 2.0, 4.0],
    vec![0.7, -0.3, 0.4, 0.5],
);

let result = couple(&sys_a, &sys_b);
println!("Mode: {:?}, Emergent: {}, Synergy: {:.4}", 
    result.mode, result.emergent, result.synergy);
```

```bash
cargo run
```

## API Reference

| Type | Description |
|------|-------------|
| `SpectralSystem` | A system with eigenvalues + dominant mode |
| `CouplingResult` | Coupled system with mode, emergence flag, synergy |
| `CouplingMode` | `Reinforcing` / `Complementary` / `Competitive` |

| Function | Returns |
|----------|---------|
| `couple(a, b)` | `CouplingResult` with emergent detection |
| `cosine_similarity(a, b)` | Alignment of two mode vectors |
| `SpectralSystem::cr()` | Conservation ratio λ₂/λₙ |
| `SpectralSystem::spectral_gap()` | Algebraic connectivity |

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
emergent-coupling = { git = "https://github.com/SuperInstance/emergent-coupling" }
```

## How It Fits

Part of the SuperInstance ecosystem:

- **[regime-detection](https://github.com/SuperInstance/regime-detection)** — Detects when coupling structure breaks
- **emergent-coupling** — Models how coupling structure emerges (this repo)

## License

MIT
