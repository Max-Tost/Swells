# 8. When groups go nonlinear

Everything in lesson 7 rests on one assumption: that the sea is a linear
superposition of independent components. That is what makes it Gaussian, which
is what makes the envelope Rayleigh, which is what gives us the whole statistical
apparatus.

At moderate steepness the assumption fails, and it fails in an interesting
direction: a group starts to pump energy into itself.

## The envelope has its own equation

We threw away the quadratic terms in lesson 1. Put the first of them back.

The technique is multiple scales. The carrier oscillates on the fast scales
$x, t$; the envelope varies on slow ones. Introduce
$\xi = \epsilon(x - c_gt)$ and $\tau = \epsilon^2t$, expand $\eta$ in powers of
$\epsilon$, and collect order by order. At $\mathcal O(\epsilon)$ you recover
linear theory. At $\mathcal O(\epsilon^2)$ you get the bound second harmonic —
the reason real waves have sharp crests and flat troughs. At
$\mathcal O(\epsilon^3)$ the solvability condition is an equation for the
envelope $A$:

$$i\left(\frac{\partial A}{\partial t} + c_g\frac{\partial A}{\partial x}\right)
+ \alpha\frac{\partial^2A}{\partial x^2} + \beta|A|^2A = 0. \tag{8.1}$$

The nonlinear Schrödinger equation. The algebra to get there is long and I am
going to skip it — Mei's *Applied Dynamics of Ocean Surface Waves* does it in
full — but the two coefficients are worth having:

$$\alpha = \frac12\frac{d^2\omega}{dk^2} = -\frac{\omega}{8k^2},
\qquad
\beta = -\frac{\omega k^2}{2}
\qquad\text{(deep water)}.$$

Read (8.1) physically. The $\alpha$ term is dispersion: it spreads a packet
out, because the components inside it travel at slightly different speeds. The
$\beta$ term is nonlinearity: a bigger local amplitude means a locally different
phase speed, which bunches the packet up. Spreading versus bunching. If they
have the right relative sign, they fight, and the fight can be won by
bunching — a group that concentrates itself.

## Benjamin and Feir

Ask what happens to a perfectly uniform wave train. Equation (8.1) has the
exact solution $A = a_0e^{i\beta a_0^2t}$ — constant amplitude, slightly shifted
frequency. That is a Stokes wave.

Perturb it with a small sideband of wavenumber $K$, linearise, and you get the
modulation frequency

$$\Omega^2 = \alpha K^2\left(\alpha K^2 - 2\beta a_0^2\right)
= \left(\frac{\omega K^2}{8k^2}\right)
\left(\frac{\omega K^2}{8k^2} - \omega k^2a_0^2\right). \tag{8.2}$$

When $\Omega^2 > 0$ the modulation just propagates. When $\Omega^2 < 0$, $\Omega$
is imaginary and the sideband **grows exponentially**. From (8.2) that happens
for

$$0 < K < 2\sqrt2\,k^2a_0,$$

with fastest growth at $K = 2k^2a_0$ and rate

$$\boxed{\ \Gamma_{\max} = \tfrac12\,\omega\,(ka_0)^2 = \tfrac12\omega\epsilon^2\ }
\tag{8.3}$$

This is the Benjamin–Feir instability, found in 1967 in a wave tank at
Feltham. A uniform train, carefully generated, spontaneously breaks up into
groups as it travels down the tank. It caused some distress at the time: it
means a perfectly regular ocean wave train is *unstable*, and the sea cannot
help but be groupy even if you try to make it otherwise.

Two numbers make it concrete.

**The scale.** At maximum growth the modulation wavelength is
$2\pi/K = \pi/(k^2a_0) = L/2\epsilon$. For $\epsilon=0.1$ that is five
wavelengths; for $\epsilon=0.02$, twenty-five. Either way, some tens of
wavelengths — which is exactly the size of a wave group. The instability picks
out the group scale on its own.

**The time.** From (8.3) the e-folding time is $2/(\omega\epsilon^2)$ seconds,
or $1/(\pi\epsilon^2)$ wave periods. At $\epsilon=0.1$ that is 32 periods,
about eight minutes for a 15 s wave. At a realistic ocean swell steepness of
$\epsilon=0.02$ it is 800 periods, three and a half hours. Slow — but the
Pacific takes a week to cross, so there is plenty of time.

## It switches off in shallow water

At finite depth the coefficients depend on $kh$, and $\beta$ **changes sign** at

$$kh = 1.363.$$

For $kh > 1.363$ the product $\alpha\beta > 0$: focusing, unstable, groups
sharpen. For $kh < 1.363$: defocusing, stable, modulations flatten out.

So a swell that is modulationally unstable in the open ocean becomes stable as
it shoals. Rogue waves are a deep-water phenomenon; the shallow ocean actively
suppresses them. For a 15 s swell, $k=0.0179$, so the switch happens at
$h = 76$ m — out on the shelf, well before the surf zone.

## Does it matter for a real sea?

The derivation above is for a single wave train with one sideband. A real sea
is a broad random field, and dispersion is constantly scrambling the phase
coherence that the instability needs. So the question is which wins: nonlinear
focusing, or dispersive spreading.

That ratio is the **Benjamin–Feir index**:

$$\mathrm{BFI} = \frac{\sqrt2\,k_p\sqrt{m_0}}{\nu}. \tag{8.4}$$

Steepness on top, bandwidth on the bottom. Above about 1, the instability
outruns the spreading, the sea develops coherent groups beyond what linear
theory predicts, and the tail of the height distribution runs *above* Rayleigh —
more very large waves than (7.2) allows.

A caution on (8.4): definitions differ. Some authors use the steepness
$k_pH_s/2$ instead of $k_p\sqrt{m_0}=k_pH_s/4$, which moves BFI by a factor of
two. So "the threshold is about 1" carries a factor-of-two ambiguity, and you
should read anyone's BFI number together with their definition.

For our default swell BFI comes out around 0.1 — nowhere near. Distant
groundswell is too flat and, after the two-swell case aside, usually too broad
for this to matter. Where it does matter is a young, steep, narrow storm sea:
short fetch, strong wind, and a very peaked spectrum. That is where the rogue
wave reports come from.

## The long wave underneath

One more nonlinear effect, and this one you can feel on any beach.

A group carries more energy than the lull between groups, and therefore more
radiation stress (lesson 9). The mean water surface responds by dipping down
beneath the big waves and rising between them. This depression is *bound* to
the group: it travels at $c_g$, not at its own free-wave speed, and it has the
period of the group — thirty to three hundred seconds.

Then the group reaches the surf zone and the short waves break. The forcing
that held the long wave in place disappears, and the long wave is released as a
free wave. It runs into the beach, reflects, and sloshes around the surf zone.

This is **surf beat**: the slow surge in and out of the water line, on a period
far longer than any individual wave. It drives rip current pulsing and it
dominates run-up on gentle beaches. In storm conditions the infragravity band
can carry more energy at the shoreline than the wind waves that made it, which
is why coastal flooding models cannot ignore it.

We do not model it here — `swells.nonlinear.bound_infragravity_amplitude` gives
only the crude scaling $\eta_{\rm ig}\sim -H_s^2/16h$ — but it is a real part of
the chain and you should know it is missing.

## What else is missing

Three honest omissions.

**Crossing seas.** Everything above is one-dimensional. Two swells crossing at
an angle behave differently: the unidirectional instability is suppressed, but
new three-dimensional patterns appear. Reports of freak waves in crossing seas
are common and the mechanism is not settled.

**Currents.** An opposing current shortens and steepens waves, which raises
$\epsilon$ and therefore the growth rate quadratically. The Agulhas has a
reputation for a reason.

**Second-order bound waves.** Real waves have sharper crests and flatter
troughs than a sine, so even without any instability the crest height
distribution runs above Rayleigh while the trough distribution runs below. This
is a systematic effect, present always, and it is not in our model.

---

*Next: the wave finds the bottom.*
