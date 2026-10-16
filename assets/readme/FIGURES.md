# MATLAB figures for the harmonic README

These three static figures were computed and rendered in MATLAB R2026b using the package at [commit `04e0574`](https://github.com/laplcebeltrami/harmonic/tree/04e0574a79c87c94c7c702a5446ab8744a7d2a57). The original package functions and the supplied simulation files were not changed. The three existing brain GIFs remain linked in the README exactly as before; they come from the separate [networkgradients repository](https://github.com/laplcebeltrami/networkgradients).

## Reproduce

From the repository root in MATLAB:

```matlab
addpath(fullfile(pwd,'examples','readme'))
generate_readme_figures
```

The [figure script](../../examples/readme/generate_readme_figures.m) writes the three PNGs to this folder and saves the numeric outputs and checks in `build/readme/`. It uses explicit figure colors and fonts, so the exported figures have a white background even when MATLAB uses a dark theme. All drawing and image export are performed in MATLAB.

## Figure 1: separate propagation and circulation

`harmonic-decomposition.png` uses the six-node network and flow strengths from the Live Script's weak-circulation setup: a four-edge cycle with weights `[2,3,3,2]`, plus two branches with weights `[5,5]`. The figure decomposes this specified flow directly; it is not an estimate from noisy time series.

`adj2simplex(C,1)` constructs its edge-only complex, `PH_boundary` supplies the boundary matrix, and `Hodge_decompose_edgeflow` computes the components. The resulting harmonic flow is `[2.5,2.5,2.5,2.5,0,0]` in the displayed edge orientation. Node positions were adjusted for readability. Every panel uses the same arrow-width scale; pale lines indicate the underlying support when a component is zero.

## Figure 2: a time-series workflow

`timeseries-to-harmonic.png` uses `simulationZ_weakcirculation.mat`, trial 1 at innovation standard deviation 1: six nodes and 200 time samples. This trial was selected before inspecting the projection.

`TS2timelagged` computes full-series lag-1 correlations. `adj_AntiSymmetric` constructs the signed directional matrix using the implementation's pairwise comparison rule. `harmonicflow` projects the resulting vector on all 15 lexicographically ordered node pairs. No known generating-edge mask or fitted threshold is applied. The node traces are standardized and vertically offset for display. The matrix colorbar shows the input flow scale; graph width encodes absolute harmonic flow and arrow direction encodes sign.

This is an example of the analysis steps, not a comparison of recovery accuracy. The `harmonicflow` entry point also avoids a name collision with the Symbolic Math Toolbox's harmonic-number function, `harmonic`.

## Figure 3: gradient, curl, and harmonic components

`hodge-components.png` uses a new deterministic six-node complex containing one filled triangle and one unfilled four-edge loop. The synthetic input is the sum of a gradient from specified node potentials, a triangular boundary flow, and a flow around the open loop. `adj2simplex(C,2)` and `PH_boundary` construct the complex; `Hodge_decompose_edgeflow` recovers the three components. The filled triangle is shaded, and all panels share one arrow-width scale.

## Checks

The executed examples verify reconstruction, agreement with the known components, orthogonality, harmonic divergence, harmonic curl, and the boundary-chain identity. Residuals for the double-precision examples were below `4e-14`. The `harmonicflow` wrapper stores its output as single precision; its matrix agreed with the double-precision decomposition within `2.8e-8` in Frobenius norm. Exact results are saved in [numerical-checks.json](numerical-checks.json).

These checks cover the README examples. They do not establish end-to-end validation of all Live Script sections, clinical validity, or superiority over the comparison methods.
