# 3. How the wind makes waves

The wind blows over flat water. Waves appear. How?

This turns out to be genuinely hard, and it stayed unsolved until 1957, when
two people solved it in the same volume of the same journal — differently, and
both correctly, because there are two mechanisms and you need both.

## The problem with the obvious answer

The obvious answer is friction: the wind drags on the surface and piles it up.
This is wrong, or at least it is answering a different question. Wind stress on
water drives a *current*, not a wave. To make a wave you need to put energy
into an oscillation, which means you need a force that pushes down on the
troughs and up on the crests — a pressure field in phase with the wave slope,
not with the wave itself.

Nothing about steady flow over a flat surface produces that. So where does it
come from?

## Phillips: resonance with the turbulence

Owen Phillips's answer (1957): the wind is turbulent, so the pressure it exerts
on the water is not smooth. It is a random, travelling field of gusts and lulls
— eddies of all sizes, sweeping downwind at roughly the wind speed.

Fourier-analyse that pressure field. It contains components of every wavenumber
$\mathbf k$, each travelling at some speed. Now: the water surface can support
a free wave of wavenumber $\mathbf k$, which travels at $c(k)=\sqrt{g/k}$. When
a pressure component happens to have the same wavenumber *and* the same speed
as a free wave, it pushes in step with that wave, cycle after cycle, and the
wave grows.

It is a driven oscillator on resonance, with the turbulence supplying the
driving. Because the forcing is random and uncorrelated with the wave it is
building, the amplitude executes a random walk, and energy — which goes as
amplitude squared — grows *linearly* in time:

$$\frac{\partial E(\mathbf k,t)}{\partial t} = \alpha_P(\mathbf k). \tag{3.1}$$

Linear growth, independent of how much wave is already there. This is the right
mechanism for getting started from a dead flat sea, and it explains the first
ripples. But linear growth is far too slow to build a storm sea, and (3.1) has
no idea that a big wave should grow faster than a small one. Something else
takes over.

## Miles: the wave builds its own forcing

John Miles's answer, in the same year, is subtler and it is the one that does
the heavy lifting.

Once there is a ripple, the airflow has to go over it. The wind profile over
the sea is roughly logarithmic, $U(z)$, increasing with height. Somewhere there
is a height $z_c$ at which the wind speed equals the wave's phase speed:

$$U(z_c) = c.$$

This is the **critical layer**. In the frame moving with the wave, the air at
$z_c$ is stationary — it sits on top of the wave, going nowhere. Above it the
air moves forwards relative to the wave, below it backwards.

That is a singular point of the linearised airflow problem, and it is where the
interesting thing happens. Miles showed that the critical layer transfers
momentum from the mean shear into the wave, and — crucially — produces a
surface pressure perturbation that is *phase shifted* relative to the surface
displacement. Not in phase with the crest, which would do no net work, but
shifted, so that the pressure is high on the upwind face and low on the
downwind face.

Now the wind is pushing the water in the direction it is already moving. Work
is done every cycle. And because the pressure perturbation is proportional to
the wave amplitude — a bigger wave distorts the airflow more, and gets pushed
harder — the growth is *exponential*:

$$\frac{\partial E(\mathbf k,t)}{\partial t} = \beta_M(\mathbf k)\,E(\mathbf k,t),
\qquad
\beta_M \sim \epsilon_a\,\omega\left(\frac{u_*}{c}\right)^2\cos^2\theta,
\tag{3.2}$$

where $\epsilon_a = \rho_{\rm air}/\rho_{\rm water} \approx 1.2\times10^{-3}$,
$u_*=\sqrt{C_d}\,U_{10}$ is the friction velocity, and $\theta$ is the angle
between the wave and the wind.

Read the physics off (3.2). The factor $\epsilon_a$ says the coupling is weak,
because air is a thousand times lighter than water — the wind has to work at
this for hours. The factor $(u_*/c)^2$ says a wave stops growing once it moves
as fast as the wind: if $c > U$, there is no critical layer, no phase shift, no
growth. **A wave cannot be driven faster than the wind that makes it.** And
$\cos^2\theta$ says waves running crosswind grow slowly and waves running
upwind not at all.

So the picture is: Phillips gets you off the ground, Miles takes over and runs
away with it, and Miles shuts off when the wave outruns the wind.

## Where the long waves come from

Except that Miles's cutoff creates a puzzle. A 25 m/s wind cannot directly
drive anything faster than 25 m/s, which by $c=gT/2\pi$ means nothing longer
than $T = 16$ s. Yet the same storm will happily produce 18 and 20 s swell.
How?

Not from the wind. From the waves themselves.

Take four wave components with wavevectors and frequencies satisfying

$$\mathbf k_1+\mathbf k_2 = \mathbf k_3+\mathbf k_4,
\qquad
\omega_1+\omega_2 = \omega_3+\omega_4. \tag{3.3}$$

These are resonance conditions, and they are exactly the conditions under which
the weak nonlinearity we discarded in lesson 1 stops being negligible. Four
waves satisfying (3.3) exchange energy with each other coherently, and over
many periods the exchange accumulates. Hasselmann worked out the transfer rate
in 1962 — a fearsome six-dimensional integral over all resonant quadruplets.

Note what (3.3) does *not* say: it says nothing about creating or destroying
energy. Four-wave interaction is conservative. It only moves energy around
within the spectrum. But where it moves it is the point: the transfer takes
energy out of the peak and dumps it on both sides, a little to higher
frequencies (where dissipation eats it) and a lot to *lower* frequencies, onto
the forward face of the spectrum.

Which means the peak walks downhill. As a sea grows, $f_p$ decreases, the waves
get longer, and eventually they get longer than the wind could ever have
driven directly. The spectrum bootstraps itself into periods the wind cannot
touch.

This is not a detail. It is the reason long-period swell exists at all, and
therefore the reason any of the rest of this project is interesting.

## The bookkeeping equation

Put the three processes together. If we track the spectrum as it evolves in
space and time, the statement is

$$\frac{\partial N}{\partial t}
+ \nabla_{\mathbf x}\cdot\left[(\mathbf c_g + \mathbf U)N\right]
+ \nabla_{\mathbf k}\cdot(\dot{\mathbf k}N)
= \frac{S_{\rm in} + S_{\rm nl} + S_{\rm ds}}{\sigma}. \tag{3.4}$$

This is the **action balance equation**, and it is what every operational wave
model on Earth integrates.

Two things about it. First, the tracked quantity is not energy but *wave
action* $N = E/\sigma$, where $\sigma = \omega - \mathbf k\cdot\mathbf U$ is
the frequency seen by someone drifting with the current. The reason is that
when waves ride on a current, energy is not conserved — the current can do work
on them — but action is. Action is to waves what particle number is to photons,
and if you have met the adiabatic invariant $E/\omega$ in mechanics, this is
the same animal.

Second, the right-hand side is the whole subject in three terms:

- $S_{\rm in}$ — wind input, the Phillips and Miles mechanisms above;
- $S_{\rm nl}$ — the four-wave transfer, which creates and destroys nothing but
  reshapes everything;
- $S_{\rm ds}$ — dissipation, overwhelmingly whitecapping. When a crest gets too
  steep it breaks, and the energy goes into turbulence. This is the term nobody
  can write down from first principles; every model uses a parameterisation
  tuned to make the answer come out right.

We are not going to solve (3.4). We are going to do something much cheaper.
Once swell leaves the storm, all three source terms essentially vanish —
$S_{\rm in}$ because there is no wind, $S_{\rm ds}$ because the waves are no
longer steep enough to break, $S_{\rm nl}$ because it is fourth order in a small
quantity. With the right-hand side zero, (3.4) is just transport, and transport
we can do by hand. That is lesson 5.

But first we need to know what spectrum the storm hands us.

---

*Next: the shape of a wind sea, and the one place where the standard results
contradict each other.*
