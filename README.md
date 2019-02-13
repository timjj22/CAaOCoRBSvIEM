# Collision-Aware and Online Compression of  Rigid Body Simulations via Integrated Error Minimization

This is a project implementation of [CAaOCoRBSvIEM](http://www.dgp.toronto.edu/projects/rigid-body-compression/).
.dmat files are easily read by both [libigl](https://github.com/libigl/libigl) and [gptoolbox](https://github.com/alecjacobson/gptoolbox) and are simple raw matrix formats.

Additional data can be provided to the compressor for contact events between objects. This allows for better collision detection and compression results.
The input data is specifed within the json file.

The **uncompressed** data format is stored in the following way:
```
 t₀  t₁  t₂  ̇⋯  tₙ

 q₀  q₀  q₀     q₀
 p₀  p₀  p₀     p₀

 q₁  q₁  q₁     q₁
 p₁  p₁  p₁     p₁

 ⋮          ⋱
 
 qₘ  qₘ  qₘ     qₘ
 pₘ  pₘ  pₘ     pₘ
```
where each of the `n` columns is a frame from the simulation, for each of the `m` objects stacked quaternion and position. Each frame fo the simulation is defined in a `(7*m + 1)` vector, stacked to include `n` frames.

The **compressed** data format is slightly different, but is stored in a similar way. Each of the `m` objects are stored in separate compressed files, and the quaternions are stored separately as well.



The rotation files: `#QcState.dmat`

Each column in this file is a keyframe for the data with the times and quaternion and angular velocity. The time in the column of the keyframe outlines the time when the keyframe ends (coefficients Xᵢ are valid for the window [tᵢ₋₁   tᵢ])
```
 t₀    t₁    ⋯    tᵢ
 α₀    α₁         αᵢ
 qx₀   qx₁        qxᵢ
 qy₀   qy₁        qyᵢ
 qz₀   qz₁        qzᵢ
```
notice how only 3 components of the quaternion are stored, they are saved in such a way that `qw = sqrt(1 - qx² - qy² - qz²)` (guaranteed to be positive).
The α is the angular velocities between the keyframes (in so-called "normalized" time, where the time at the start and end of the keyframes are 0 and 1 respectively), this way multiple spins around an axis can be encoded. In order to avoid numerical issues with the quaternion packages, simple slerps are used if the angular velocity term is too small.



The position files: `#PcState.dmat`

Each column in the file is a keyframe with the times and positions and velocities. The positions are reconstructed from the following quadratic: `p(t) = αt² + ̱βt + γ`, again in normalized time. The data is stored as shown:
```
 tᵢ    ⋯    tᵢ₊₁
 βᵢ         βᵢ₊₁
 γᵢ₊₁       γᵢ₊₂
```

Endpoints at the end of the current keyframe and the start of the next are the same, and so the quadratic can be reconstructed as `αᵢ₊₁ = γᵢ₊₂ - γᵢ₊₁ - βᵢ₊₁`
