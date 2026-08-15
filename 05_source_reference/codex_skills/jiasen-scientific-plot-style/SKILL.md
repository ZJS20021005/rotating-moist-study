---
name: jiasen-scientific-plot-style
description: Use for scientific figures in this workspace when the user asks to画图, restyle plots, keep a fixed paper-like format, or says to follow the previously agreed plot style. Applies white background, thick black frame, inward ticks, STIX fonts, fixed case colors, and a fixed single-panel aspect ratio that must also be preserved inside multi-panel figures.
---

# Jiasen Scientific Plot Style

Use this skill for all scientific plots unless the user explicitly asks for a different style.

## Core Style

- White background.
- Final figures should use the approved `q_{\rm rms}(z)` profile standard by default: clean white background, moderately thick black frame, axes-box aspect ratio preserved, inward ticks, and saturated colors. Use a different example style only when the user explicitly says to follow another provided figure.
- Ticks point inward on all four sides.
- Turn on minor ticks.
- Use Times New Roman fonts:
  - `font.family = Times New Roman`
  - keep math text in a serif style that matches as closely as possible
- Use thicker lines than matplotlib defaults.
- Prefer clear legends that do not overlap important curves.

Recommended baseline rcParams:

```python
mpl.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "axes.linewidth": 4.5,
    "axes.labelsize": 24,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 10,
    "figure.dpi": 180,
    "savefig.dpi": 300,
})
```

Recommended axis styling:

```python
for spine in ax.spines.values():
    spine.set_linewidth(4.5)
ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
ax.minorticks_on()
ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
```

Treat the frame and font values above as the approved visual standard, not
loose suggestions: frame `4.5`, axis-title size `24`, tick-label size `13`,
and legend size `10`. At 300 dpi with the fixed `6.5 x 5.84` reference canvas
and axes rectangle, the frame is about 19 pixels thick, matching the approved
reference after display scaling. Use the same values for figures with and
without legends.

## Fixed Panel Aspect Ratio

This is mandatory.

- Standard single-panel outer canvas: `6.5 x 5.84` inches. The black axes box itself must keep the `6.5:5.2` width/height ratio.
- The required width/height ratio refers to the **black axes frame (the plotting box)**, not the outer PNG canvas.
- The width/height ratio of every panel's **axes box** in a multi-panel figure must match the single-panel ratio.
- Do not squash or stretch panels just because the figure is a grid.

Practical rule:

- preserve the axes-box ratio with something like:

```python
ax.set_box_aspect(5.2 / 6.5)
```

- then choose the outer figure size large enough to hold legends, labels, and colorbars without changing the axes-box ratio.

Use these defaults:

- `1x1`: `6.5 x 5.84`
- `1x2`: `13.0 x 5.84`
- `2x2`: `13.0 x 11.68`
- `4x3`: `19.5 x 23.36`

Small adjustments for margins or colorbars are allowed, but keep each panel close to the same ratio.

## Colors and Line Conventions

Keep case colors consistent across related figures.

Preferred palette used repeatedly in this workspace comes from:

`E:\moist RB\post\vprofile`

Concrete colors:

```python
[
    (0.74, 0.14, 0.18),
    (0.93, 0.32, 0.23),
    (0.96, 0.58, 0.19),
    (0.38, 0.75, 0.91),
    (0.16, 0.56, 0.80),
    (0.11, 0.44, 0.71),
]
```

Rules:

- Keep the same case-color mapping across all figures in a set.
- Use black dashed lines for `dry` when it appears as a reference/control.
- Do not invent a new palette if an existing case color mapping already exists in the working set.
- For current rotating-moist Ra comparison figures, use the high-saturation case colors unless the user asks otherwise:

```python
RA_COLORS = {
    1.0e6: (0.00, 0.65, 0.12),   # saturated green
    5.0e6: (0.58, 0.05, 0.90),   # saturated purple
    5.0e7: (1.00, 0.00, 0.00),   # saturated red
    1.0e8: (0.00, 0.00, 1.00),   # saturated blue
}
```

Use filled markers by default for these figures, including the `Ra=5e7` red points, unless the user explicitly asks for hollow markers.

## Legends

- Do not let legends cover important lines.
- Prefer a single clean legend over multiple overlapping legends.
- Legends should be frameless by default.
- Only include scalar summary values in the legend if the user explicitly asks.
- If legends become crowded, move them to an empty region inside the axes before considering a figure layout change.

## Labels

- Use math notation for scientific quantities.
- Keep axis labels concise.
- Do not write case labels as literal `q_bot`; render the subscript in math notation, for example `r"$q_{\mathrm{bot}}=0.5q_s$"`.
- For wall Nusselt labels, use exactly:
  - `Nu_low`
  - `Nu_top`
- For volume-averaged Nusselt based on `avgvar` last column, use:

```python
r"$Nu_{\mathrm{vol}}=\left\langle \sqrt{RaPr}\,wT-\partial_z T \\right\rangle_V$"
```

Use shorter labels on axes when needed, and move long expressions to captions or summary text.

## Scientific notation

- For a vertical axis whose characteristic magnitude is very small or very
  large, use a shared scientific multiplier instead of repeating powers of
  ten in every tick label. The project default is Matplotlib power limits
  `(-2, 2)` with STIX math text:

```python
ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(13)
```

- Put a left-y-axis multiplier at the upper left of the frame. For a dual-y
  plot, put the right-y-axis multiplier at the upper right. Keep both axes,
  tick marks, tick labels, and axis titles black; use curve colors plus a
  frameless legend to distinguish quantities.
- If an x axis requires scientific notation, place its multiplier at the
  lower right. Do not repeat `10^n` on every tick.

## Time Series Plots

- Use the same thick-frame style.
- Single-panel time series can use `7.0 x 5.4` when a little extra width helps readability.
- For comparison plots, keep all cases on the same axes unless the user explicitly asks to split them.
- If plotting bottom/top wall quantities together, use a `1x2` figure and preserve panel ratio.

## Default Profile Plot Standard

This is now the default style for profile figures and, unless the user explicitly says otherwise, for future scientific plots in this project. It is based on the previously approved `q_{\rm rms}(z)` profile figure.

Use these exact defaults when drawing vertical profiles such as `q_{\rm rms}(z)`, `w_{\rm rms}(z)`, `u_{\rm rms}(z)`, `m(z)`, `RH(z)`, and related profile comparisons:

```python
mpl.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "axes.linewidth": 4.5,
    "axes.labelsize": 24,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 10,
    "figure.dpi": 180,
    "savefig.dpi": 300,
})

fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])

for spine in ax.spines.values():
    spine.set_linewidth(4.5)
ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
ax.minorticks_on()
ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
ax.set_box_aspect(5.2 / 6.5)
```

For multi-case Ek profile plots, use:

```python
colors = plt.cm.turbo(np.linspace(0.06, 0.94, len(eks)))
ax.plot(profile, z, color=color, lw=3.5, label=ek_label(ek))
```

Preferred legend for profile comparison figures:

```python
ax.legend(
    title=r"$Ek$",
    frameon=False,
    ncol=2,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
)
fig.subplots_adjust(left=0.14, right=0.72, bottom=0.14, top=0.96)
```

If the user asks for "不要 label / 单独 label 图", remove the legend from the main figure and create a standalone legend figure using the same line colors and line widths.

When producing both a with-legend and a no-legend version of the same figure, keep the plotted axes box physically identical. Do not use different `subplots_adjust` margins that change the axes width or height. Use an explicit shared axes rectangle for both versions, for example:

```python
FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]

fig = plt.figure(figsize=FIGSIZE, facecolor="white")
ax = fig.add_axes(AX_RECT)
```

Put the legend in an empty region inside the axes, or create a standalone legend figure. Do not shrink the axes to make room for a legend. Save paired with/without-legend figures without `tight_layout` and without `bbox_inches="tight"`. This keeps the canvas and all four frame coordinates identical.

Only deviate from this profile standard when:

- the user explicitly says not to use this style;
- the user provides another example figure and asks to match that example;
- the plot type genuinely requires a different layout, such as a joint PDF with colorbar or an x-z field with streamlines.

## Profile Plots

- Default vertical profile outer canvas: `6.5 x 5.84`; preserve a `6.5:5.2` black axes box.
- Use consistent line widths around `3.5` for normal profile curves.
- For paired speed plots, keep `u_h` solid and `w` dashed unless the user specifies another convention.
- For profile comparison figures like `u_{\rm rms}(z)` or `q_{\rm rms}(z)`, use the approved high-saturation `turbo` profile style unless the user asks otherwise:

```python
colors = plt.cm.turbo(np.linspace(0.06, 0.94, n_cases))
```

  - Thick black frame: `axes.linewidth = 4.5`.
  - Long inward major ticks on all sides: `length=12`, `width=1.2`, `top=True`, `right=True`.
  - Minor ticks on all sides: `length=6`, `width=1.0`.
  - Profile line width about `3.5` for final figures, at least `3.0` for crowded multi-case figures.
  - Preserve the fixed panel aspect with `ax.set_box_aspect(5.2 / 6.5)`.
  - Use Times New Roman/STIX math labels.
  - Use saturated lines without hollow markers for profile curves.
- Add reference lines only when requested or when they are a stable part of the interpretation, such as:
  - `RH = 0.98`
  - `kurtosis = 3`

## Joint PDF Plots

- Preserve panel aspect ratio.
- Use the agreed yellow-magenta-blue palette when matching the previously shown PDF style.
- Keep the same global axis range and color scale across all cases in the same comparison figure.
- For multi-case, multi-height PDFs:
  - rows = cases
  - columns = heights
- Add dashed zero lines when they help interpretation.
- Prefer a single shared colorbar.

## Output and Workflow Rules

- Save figures directly to the user-specified output directory.
- If the user has already fixed an output directory for a figure set, continue using it.
- When updating a figure, overwrite the old file instead of creating unnecessary variants unless the user asks for both.
- If a previous plot violated the style, restyle it instead of arguing about the old version.

## Historical Base Style

This skill should follow the historical plotting style used in:

`E:\moist RB\post\vprofile`

Keep these baseline traits from that directory:

- single-panel outer canvas `6.5 x 5.84`, with a `6.5:5.2` black axes box
- white background
- inward ticks
- legend without frame
- warm-to-cool fixed palette across cases

Then apply the newer stricter rules from this skill on top of that base:

- thicker black frame than the old scripts
- consistent case naming and legend rules
- fixed panel ratio inside multi-panel figures
- Times New Roman font family

## When Unsure

Default to this style.

Only deviate if the user explicitly asks for:

- a different journal style
- a lightweight diagnostic plot
- a screenshot-like quick check rather than a publication-style figure
