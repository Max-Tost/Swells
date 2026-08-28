"""Wave groups: envelopes, height statistics, and how many waves are in a set.

Derived in lessons/07-why-waves-come-in-sets.md.

This is the module that answers the question the project was built around. Two
independent routes to the same number are provided, and they are meant to be
compared:

  * theory  -- the bivariate Rayleigh model of Kimura (1980), which predicts
               the mean run length from the spectrum alone;
  * counting -- take a synthesised record, find the waves, count the runs.

If they disagree by more than the sampling error, one of them is wrong.
"""

import numpy as np
from scipy.signal import hilbert
from scipy.special import i0e

from .spectra import mean_frequency, mean_period, moment
from .util import trapz


# ---------------------------------------------------------------------------
# Envelope
# ---------------------------------------------------------------------------

def envelope(eta):
    """The slowly varying amplitude R(t) via the Hilbert transform.

    Write eta(t) = R(t) cos(omega_bar t + theta(t)). The analytic signal
    eta + i H[eta] has modulus R. For a narrow-band sea this is exactly the
    curve you would draw through the crests -- and a set is a hump in it.
    """
    return np.abs(hilbert(np.asarray(eta, dtype=np.float64)))


def envelope_correlation(f, S, tau=None):
    """kappa: the correlation of the complex envelope at lag tau.

        rho(tau) = |int S(f) exp(-2 pi i (f - f_bar) tau) df| / m_0

    Evaluated by default at tau = the mean wave period, because that is the lag
    between one wave and the next -- which is what governs whether a big wave
    is followed by another big wave.

    Note the modulus of the *complex* integral. Taking only the cosine part, as
    is sometimes written, throws away the sine term and underestimates kappa
    for asymmetric spectra.
    """
    f = np.asarray(f, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)
    m0 = moment(f, S, 0)
    if m0 <= 0:
        return 0.0
    if tau is None:
        tau = mean_period(f, S)
    f_bar = mean_frequency(f, S)
    phase = -2.0j * np.pi * (f - f_bar) * tau
    rho = trapz(S * np.exp(phase), f) / m0
    return float(np.clip(np.abs(rho), 0.0, 0.999999))


# ---------------------------------------------------------------------------
# Height statistics
# ---------------------------------------------------------------------------

def rayleigh_exceedance(H, Hm0):
    """P(height > H) = exp(-2 (H/H_m0)^2), for a narrow-band Gaussian sea."""
    return np.exp(-2.0 * (np.asarray(H, dtype=np.float64) / Hm0) ** 2)


def height_statistics(Hm0):
    """The Rayleigh family of height statistics, all fixed by H_m0 alone.

    Ratios to H_m0: mean 0.626, rms 0.707, H_1/3 1.000, H_1/10 1.272,
    H_1/100 1.668. That last one is the useful one -- the biggest wave in a
    couple of hours is roughly 1.7 times the significant height, which is why
    a "two metre" day has three-and-a-half metre sets in it.
    """
    return {
        "Hm0": Hm0,
        "H_mean": 0.6267 * Hm0,
        "H_rms": Hm0 / np.sqrt(2.0),
        "H_third": Hm0,
        "H_tenth": 1.2716 * Hm0,
        "H_hundredth": 1.6684 * Hm0,
    }


# ---------------------------------------------------------------------------
# Run lengths: the bivariate Rayleigh model
# ---------------------------------------------------------------------------

def bivariate_rayleigh_pdf(x1, x2, kappa):
    """Joint pdf of two successive normalised heights x = H / H_rms.

        p = (4 x1 x2 / (1-k^2)) I_0(2 k x1 x2/(1-k^2))
                              exp[-(x1^2 + x2^2 - 2 k x1 x2)/(1-k^2)]

    Note the exponent has been folded together with the I_0 argument so that
    nothing overflows: we use the exponentially scaled i0e, and the leftover
    exponent -(x1-x2)^2-ish term stays negative.

    At kappa = 0 this factorises into two independent Rayleighs, as it must.
    """
    d = 1.0 - kappa**2
    z = 2.0 * kappa * x1 * x2 / d
    return (4.0 * x1 * x2 / d) * i0e(z) * np.exp(-(x1**2 + x2**2 - 2.0 * kappa * x1 * x2) / d)


def p22(kappa, threshold_ratio=np.sqrt(2.0), n=400, span=6.0):
    """P(next wave also exceeds the threshold | this one did).

    threshold_ratio is the cut expressed in units of H_rms. The default
    sqrt(2) is H_c = H_m0, the usual definition of a "big" wave, since
    H_m0 = sqrt(2) H_rms.

        p22 = P(x1 > xc and x2 > xc) / P(x1 > xc),  P(x > xc) = exp(-xc^2)

    Evaluated by 2-D quadrature. Analytic only in the kappa = 0 case, where it
    collapses to exp(-xc^2) -- successive waves independent, no grouping.
    """
    xc = float(threshold_ratio)
    kappa = float(np.clip(kappa, 0.0, 0.999999))

    x = np.linspace(xc, xc + span, n)
    X1, X2 = np.meshgrid(x, x, indexing="ij")
    joint = trapz(trapz(bivariate_rayleigh_pdf(X1, X2, kappa), x, axis=1), x)
    marginal = np.exp(-(xc**2))
    return float(np.clip(joint / marginal, 0.0, 0.999999))


def run_length_distribution(p, j_max=30):
    """P(run of exactly j waves) = p^(j-1) (1-p). Geometric.

    Geometric because the bivariate model is Markov: whether wave j+1 is big
    depends on wave j and nothing earlier.
    """
    j = np.arange(1, j_max + 1)
    return j, p ** (j - 1) * (1.0 - p)


def mean_run_length(p):
    """E[j] = 1/(1-p). Waves per set."""
    return np.inf if p >= 1.0 else 1.0 / (1.0 - p)


def set_interval(mean_j, T_mean, threshold_ratio=np.sqrt(2.0)):
    """Mean time from the start of one set to the start of the next.

    A cycle is one run of big waves plus the lull that follows. In a stationary
    sea the fraction of waves that are big is q = exp(-xc^2), and the big ones
    all live inside the runs, so

        E[cycle in waves] = E[run] / q,   set interval = E[run] T_mean / q.

    With the usual threshold H_c = H_m0 we have q = e^-2 = 0.135, so a cycle is
    about 7.4 runs long. Two-wave sets at a 12 s period then come round every
    three minutes, which is what you actually sit through.
    """
    q = np.exp(-(float(threshold_ratio) ** 2))
    return mean_j * T_mean / q


def group_statistics(f, S, threshold_ratio=np.sqrt(2.0)):
    """Everything about grouping, predicted from the spectrum alone."""
    kappa = envelope_correlation(f, S)
    p = p22(kappa, threshold_ratio)
    mean_j = mean_run_length(p)
    T_mean = mean_period(f, S)
    return {
        "kappa": kappa,
        "p22": p,
        "mean_run_length": mean_j,
        "sd_run_length": np.sqrt(p) / (1.0 - p) if p < 1 else np.inf,
        "mean_period": T_mean,
        "set_interval": set_interval(mean_j, T_mean, threshold_ratio),
    }


# ---------------------------------------------------------------------------
# Counting, for comparison
# ---------------------------------------------------------------------------

def zero_crossing_waves(t, eta):
    """Split a record into individual waves at zero up-crossings.

    A wave runs from one up-crossing of the mean level to the next, and its
    height is the largest crest minus the deepest trough in between. This is
    the definition a buoy report uses.

    Returns (heights, periods).
    """
    eta = np.asarray(eta, dtype=np.float64) - np.mean(eta)
    t = np.asarray(t, dtype=np.float64)

    up = np.where((eta[:-1] <= 0) & (eta[1:] > 0))[0]
    if up.size < 2:
        return np.array([]), np.array([])

    # Refine each crossing by linear interpolation.
    frac = -eta[up] / (eta[up + 1] - eta[up])
    t_cross = t[up] + frac * (t[up + 1] - t[up])

    heights, periods = [], []
    for i in range(up.size - 1):
        seg = eta[up[i]:up[i + 1] + 1]
        if seg.size < 2:
            continue
        heights.append(seg.max() - seg.min())
        periods.append(t_cross[i + 1] - t_cross[i])
    return np.array(heights), np.array(periods)


def count_runs(heights, threshold):
    """Lengths of the runs of consecutive waves exceeding the threshold."""
    big = np.asarray(heights) > threshold
    runs, current = [], 0
    for b in big:
        if b:
            current += 1
        elif current:
            runs.append(current)
            current = 0
    if current:
        runs.append(current)
    return np.array(runs)


def measured_group_statistics(t, eta, threshold=None):
    """Grouping measured by counting, for comparison with the theory above."""
    H, T = zero_crossing_waves(t, eta)
    if H.size == 0:
        return {"n_waves": 0}
    Hm0 = 4.0 * np.std(eta)
    thr = Hm0 if threshold is None else threshold
    runs = count_runs(H, thr)
    return {
        "n_waves": H.size,
        "Hm0_from_record": Hm0,
        "H_mean": float(np.mean(H)),
        "H_tenth": float(np.mean(np.sort(H)[-max(1, H.size // 10):])),
        "mean_period": float(np.mean(T)),
        "n_groups": runs.size,
        "mean_run_length": float(np.mean(runs)) if runs.size else 0.0,
        "runs": runs,
    }


def beat_period(f1, f2):
    """Two frequencies beat with envelope period 1/|f1 - f2|.

    The simplest possible set: 14 s and 15 s swell together give a beat every
    1/(1/14 - 1/15) = 210 s, about three and a half minutes, which is exactly
    the set interval surfers quote.
    """
    return 1.0 / abs(f1 - f2)
