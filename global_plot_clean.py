
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

# --------------------
# Pfade & Parameter
# --------------------
plot_dir = "/home/wernechi/NCT/results/plots_clean"
os.makedirs(plot_dir, exist_ok=True)

frequency_bands = [
    "alpha", "betalow", "betahigh", "delta",
    "gammalow", "gammamid", "gammahigh", "theta"
]

pathtoglob = "/home/wernechi/NCT/results/01_NCT-final/glob_csv/"
pathtomeandiff = "/home/wernechi/NCT/results/02_NCT-tested/00_meandiff/mean_diff_csv/"

sns.set_context("talk", font_scale=1.4)

# --------------------
# vordefinierte p-Werte (aus ANOVA)
# --------------------
pvals_stab = {
    "alpha": 0.0001, "betalow": 0.002, "betahigh": 0.208, "delta": 0.002,
    "gammalow": 0.390, "gammamid": 0.001, "gammahigh": 0.121, "theta": 0.0001
}
pvals_energy = {
    "alpha": 0.0001, "betalow": 0.0001, "betahigh": 0.0001, "delta": 0.009,
    "gammalow": 0.440, "gammamid": 0.0001, "gammahigh": 0.0001, "theta": 0.015
}

def star_label(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return None

# --------------------
# 4er-Blöcke bilden
# --------------------
quads_simple = [frequency_bands[i:i+4] for i in range(0, len(frequency_bands), 4)]

# =====================
# STABILITY – 4 nebeneinander
# =====================
with open(os.path.join(plot_dir, "stability_medians_quartiles.txt"), "w") as f_out:
    for qi, quad in enumerate(quads_simple, start=1):
        fig, axes = plt.subplots(1, 4, figsize=(24, 6), sharey=False)
        for j in range(len(quad), 4):
            fig.delaxes(axes[j])  # falls weniger als 4 Bänder im letzten Block

        for ax, band in zip(axes, quad):
            fcsv = f"{pathtoglob}{band}_wm_0B_2B_glob_T1.csv"
            df_stab = pd.read_csv(fcsv)
            df_stab["Subject"] = range(1, len(df_stab) + 1)

            # Mediane & Quartile berechnen und speichern
            for cond in ["Stab A", "Stab B"]:
                median = df_stab[cond].median()
                q1 = df_stab[cond].quantile(0.25)
                q3 = df_stab[cond].quantile(0.75)
                f_out.write(f"{band} - {cond}: Median={median:.3f}, Q1={q1:.3f}, Q3={q3:.3f}\n")

            long_df = pd.melt(
                df_stab, id_vars="Subject",
                value_vars=["Stab A", "Stab B"],
                var_name="Condition", value_name="Stability"
            )
            order_stab = ["Stab A", "Stab B"]
            long_df["Condition"] = pd.Categorical(long_df["Condition"], categories=order_stab, ordered=True)

            # Violinplot + Punkte
            sns.violinplot(
                x="Condition", y="Stability", data=long_df, ax=ax,
                inner="quartile", color="#afdebd", order=order_stab, cut=0
            )
            sns.stripplot(
                x="Condition", y="Stability", data=long_df, ax=ax,
                color="#080808", size=2, jitter=True, alpha=0.6, order=order_stab
            )

            ax.set_title(band, pad=40, fontsize=20)
            ax.set_xlabel("Condition", fontsize=22)
            ax.set_ylabel("Stability", fontsize=22)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["0-back", "2-back"], fontsize=18)

            # Signifikanzsternchen
            p = pvals_stab.get(band, None)
            if p is not None:
                label = star_label(p)
                if label:
                    pairs = [("Stab A", "Stab B")]
                    annot = Annotator(ax, pairs, data=long_df, x="Condition", y="Stability", order=order_stab)
                    annot.configure(test=None, text_format="simple", loc="outside",
                                line_offset=0.20, fontsize=18, verbose=0)
                    annot.set_custom_annotations([label])
                    annot.annotate()

        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"stability_violinplots_quad{qi}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

# =====================
# ENERGY – 4 nebeneinander
# =====================
with open(os.path.join(plot_dir, "energy_medians_quartiles.txt"), "w") as f_out:
    for qi, quad in enumerate(quads_simple, start=1):
        fig, axes = plt.subplots(1, 4, figsize=(24, 6), sharey=False)
        for j in range(len(quad), 4):
            fig.delaxes(axes[j])

        for ax, band in zip(axes, quad):
            fcsv = f"{pathtoglob}{band}_wm_0B_2B_glob_T1.csv"
            df_en = pd.read_csv(fcsv)
            df_en["Subject"] = range(1, len(df_en) + 1)

            # Mediane & Quartile berechnen und speichern
            for cond in ["Ener A-B", "Ener B-A"]:
                median = df_en[cond].median()
                q1 = df_en[cond].quantile(0.25)
                q3 = df_en[cond].quantile(0.75)
                f_out.write(f"{band} - {cond}: Median={median:.3f}, Q1={q1:.3f}, Q3={q3:.3f}\n")

            long_df = pd.melt(
                df_en, id_vars="Subject",
                value_vars=["Ener A-B", "Ener B-A"],
                var_name="Condition", value_name="Energy"
            )
            order_ener = ["Ener A-B", "Ener B-A"]
            long_df["Condition"] = pd.Categorical(long_df["Condition"], categories=order_ener, ordered=True)

            sns.violinplot(
                x="Condition", y="Energy", data=long_df, ax=ax,
                inner="quartile", color="#dba4da", order=order_ener, cut=0
            )
            sns.stripplot(
                x="Condition", y="Energy", data=long_df, ax=ax,
                color="#080808", size=2, jitter=True, alpha=0.6, order=order_ener
            )

            ax.set_title(band, pad=40, fontsize=20)
            ax.set_xlabel("Condition", fontsize=22)
            ax.set_ylabel("Transition Energy", fontsize=22)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["2-back→0-back", "0-back→2-back"], fontsize=18)

            # Signifikanzsternchen
            p = pvals_energy.get(band, None)
            if p is not None:
                label = star_label(p)
                if label:
                    pairs = [("Ener A-B", "Ener B-A")]
                    annot = Annotator(ax, pairs, data=long_df, x="Condition", y="Energy", order=order_ener)
                    annot.configure(test=None, text_format="simple", loc="outside",
                                line_offset=0.20, fontsize=18, verbose=0)
                    annot.set_custom_annotations([label])
                    annot.annotate()

        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"energy_violinplots_quad{qi}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

