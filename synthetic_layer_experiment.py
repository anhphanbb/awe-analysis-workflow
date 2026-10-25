"""Standalone synthetic two-layer experiment. Requires numpy, scipy, matplotlib.

Run: python synthetic_layer_experiment.py --iterations 256 --output results_256
No AWE data, network access, or other project scripts are required.
"""
from pathlib import Path
import argparse
import json
import warnings
import numpy as np
from scipy.ndimage import map_coordinates
import matplotlib.pyplot as plt


def run_experiment(iterations=256, output='synthetic_results', noise_sigma=0.06,
                   seed=42, report_every=8):
    """Run a deterministic experiment; save metrics, snapshots and plots.

    Known velocities are supplied, not estimated. Each iteration updates fast,
    then slow, using all 50 frames. Unsupported updates retain previous values.
    Returns a dictionary so notebook users can inspect arrays and metrics.
    """
    if iterations < 1 or noise_sigma < 0 or report_every < 1:
        raise ValueError('iterations/report_every must be positive; noise_sigma >= 0')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    y, x = np.indices((160, 160), dtype=float)
    times = (np.arange(50) - 25) * 1.1
    va, vc = np.array([.12, .07]), np.array([.48, -.40])
    roi = np.s_[40:120, 40:120]

    def air(x, y):
        return 8 + .65*np.sin(2*np.pi*(x+.40*y)/23) + .30*np.cos(2*np.pi*(y-.25*x)/39)

    def cloud(x, y):
        return (2.5*np.exp(-((x-73)**2/100+(y-88)**2/150))
                + 1.8*np.exp(-((x-94)**2/65+(y-68)**2/100))
                + 1.3*np.exp(-((x-60)**2/55+(y-65)**2/55))
                + .9*np.exp(-((x-82)**2/600+(y-80)**2/420)))

    def warp(a, v, dt):
        coords = [y-v[1]*dt, x-v[0]*dt]
        weight = map_coordinates(np.isfinite(a).astype(float), coords, order=1,
                                 mode='constant', cval=0)
        value = map_coordinates(np.nan_to_num(a), coords, order=1,
                                mode='constant', cval=0)
        return np.where(weight > .999, value/np.maximum(weight, 1e-12), np.nan)

    def predict(template, velocity):
        return np.array([warp(template, velocity, t) for t in times])

    def stack(residual, velocity):
        aligned = np.array([warp(a, velocity, -t) for a, t in zip(residual, times)])
        support = np.isfinite(aligned).sum(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            template = np.nanmedian(aligned, axis=0)
        template[support < 10] = np.nan
        return template, support

    true_slow, true_fast = air(x, y), cloud(x, y)
    rng = np.random.default_rng(seed)
    observed = (np.array([air(x-va[0]*t, y-va[1]*t) for t in times])
                + np.array([cloud(x-vc[0]*t, y-vc[1]*t) for t in times])
                + rng.normal(0, noise_sigma, (50, 160, 160)))
    slow, slow_support = stack(observed, va)
    fast = np.zeros_like(slow)
    fast_support = np.zeros_like(slow_support)
    save_at = {0, iterations} | {2**i for i in range(iterations.bit_length()) if 2**i <= iterations}
    snapshots, metrics = {}, []

    def record(n, slowseq, fastseq):
        # Remove just one shared additive constant for truth-error evaluation.
        # Estimated templates themselves are not adjusted by this correction.
        offset = np.mean((slow-true_slow)[roi])
        ea, ec = (slow-offset-true_slow)[roi], (fast+offset-true_fast)[roi]
        residual = (observed-slowseq-fastseq)[:, 40:120, 40:120]
        if not (np.isfinite(ea).all() and np.isfinite(ec).all() and np.isfinite(residual).all()):
            raise RuntimeError('Fixed evaluation region lost valid coverage')
        m = dict(iteration=n, slow_contrast_rmse=float(np.sqrt(np.mean(ea**2))),
                 fast_contrast_rmse=float(np.sqrt(np.mean(ec**2))),
                 observation_rmse=float(np.sqrt(np.mean(residual**2))),
                 shared_offset_for_truth_comparison=float(offset),
                 slow_updated_pixels=int(np.count_nonzero(slow_support >= 10)),
                 fast_updated_pixels=int(np.count_nonzero(fast_support >= 10)))
        metrics.append(m)
        if n in save_at:
            snapshots[n] = (slow.copy(), fast.copy())
        # Save progress so a long run's completed metrics remain available.
        if n % report_every == 0 or n == iterations:
            (out/'metrics.json').write_text(json.dumps(metrics, indent=2), encoding='utf-8')
            print(f"Iteration {n:4d}: slow={m['slow_contrast_rmse']:.5f} "
                  f"fast={m['fast_contrast_rmse']:.5f} observation={m['observation_rmse']:.5f}", flush=True)

    slowseq = predict(slow, va)
    record(0, slowseq, predict(fast, vc))
    for n in range(1, iterations+1):
        candidate, fast_support = stack(observed-slowseq, vc)
        fast = np.where(np.isfinite(candidate), candidate, fast)
        offset = np.nanmedian(fast)
        fast -= offset
        slow += offset
        fastseq = predict(fast, vc)
        candidate, slow_support = stack(observed-fastseq, va)
        slow = np.where(np.isfinite(candidate), candidate, slow)
        slowseq = predict(slow, va)
        record(n, slowseq, fastseq)

    np.savez_compressed(out/'experiment.npz', observed=observed, true_slow=true_slow,
                        true_fast=true_fast, times=times, slow_velocity=va,
                        fast_velocity=vc, iterations=list(snapshots),
                        estimated_slow=np.array([a for a,c in snapshots.values()]),
                        estimated_fast=np.array([c for a,c in snapshots.values()]))
    (out/'config.json').write_text(json.dumps(dict(iterations=iterations, seed=seed,
         noise_sigma=noise_sigma, frames=50, cadence_seconds=1.1, min_support=10,
         evaluation_region='[40:120,40:120]', boundary_rule='retain previous unsupported values'), indent=2))
    fig, ax = plt.subplots(figsize=(9,5), layout='constrained')
    for key, label in [('slow_contrast_rmse','Slow contrast error'),
                       ('fast_contrast_rmse','Fast contrast error'),
                       ('observation_rmse','Observation reconstruction error')]:
        ax.plot([m['iteration'] for m in metrics], [m[key] for m in metrics], label=label)
    ax.axhline(noise_sigma, color='gray', ls='--', label='Input noise sigma (reference)')
    ax.set(xlabel='Completed iterations', ylabel='RMSE (synthetic units)',
           title='Fixed central region; same observations at every iteration')
    ax.legend(); ax.grid(alpha=.25)
    fig.savefig(out/'errors.png', dpi=150); plt.close(fig)
    for n, (a,c) in snapshots.items():
        offset = np.mean((a-true_slow)[roi])
        fig, axes = plt.subplots(2,3, figsize=(12,8), layout='constrained')
        arrays = [true_slow,a,a-offset-true_slow,true_fast,c,c+offset-true_fast]
        labels = ['True slow','Estimated slow','Slow contrast error','True fast','Estimated fast','Fast contrast error']
        for k,(ax,arr,label) in enumerate(zip(axes.flat,arrays,labels)):
            cmap,lo,hi = ('RdBu_r',-1.2,1.2) if k%3==2 else (('gray',6.9,9.1) if k<3 else ('magma',-.4,4))
            im=ax.imshow(arr[roi],origin='lower',extent=[40,120,40,120],cmap=cmap,vmin=lo,vmax=hi)
            ax.set(title=label,xlabel='x (pixels)',ylabel='y (pixels)')
            fig.colorbar(im,ax=ax,shrink=.75)
        fig.suptitle(f'Iteration {n} | Fixed scales | Truth errors remove one shared constant')
        fig.savefig(out/f'iteration_{n:04d}.png',dpi=120);plt.close(fig)
    print(f'Results saved to {out.resolve()}')
    return dict(metrics=metrics, snapshots=snapshots, observed=observed,
                true_slow=true_slow, true_fast=true_fast, output=out)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--iterations', type=int, default=256)
    parser.add_argument('--output', default='synthetic_results_256')
    parser.add_argument('--noise-sigma',type=float,default=.06)
    parser.add_argument('--seed',type=int,default=42)
    args=parser.parse_args()
    run_experiment(args.iterations,args.output,args.noise_sigma,args.seed)
