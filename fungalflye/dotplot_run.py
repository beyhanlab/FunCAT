import os
import matplotlib.pyplot as plt
from pathlib import Path

from .dotplot.FastaGenome import FastaGenome
from .dotplot.Genome import MemGeneSet, spanningLocus
from .dotplot.GenomeCoord import GenomeCoord
from .dotplot.Locus import Locus
from .dotplot.MUMmerTools import NucmerMap, ClickCoords


def run_dotplot(genome1_path, genome2_path, output_dir):

    genome1_path = Path(genome1_path)
    genome2_path = Path(genome2_path)

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    genome1_name = genome1_path.stem
    genome2_name = genome2_path.stem

    print(f"\n[funcat] Comparing {genome1_name} vs {genome2_name}")

    alignment_output = output_dir / f"{genome1_name}_{genome2_name}"

    delta_file = f"{alignment_output}.delta"
    coords_file = f"{alignment_output}.coords"

    os.system(
        f"nucmer -p {alignment_output} {genome1_path} {genome2_path}"
    )

    os.system(
        f"show-coords -r -c -l -d {delta_file} > {coords_file}"
    )

    genome1 = FastaGenome(open(genome1_path, "rt"), genome1_name)
    genome2 = FastaGenome(open(genome2_path, "rt"), genome2_name)

    comp = NucmerMap.from_coords(
        coords_file,
        genome1,
        genome2,
        GenomeCoord(genome1),
        GenomeCoord(genome2)
    )

    fig = comp.plot()

    click_coords = ClickCoords(comp)
    fig.canvas.mpl_connect('button_press_event', click_coords)

    plt.title(f"{genome1_name} vs {genome2_name}")

    output_plot = output_dir / f"{genome1_name}_{genome2_name}_dotplot.pdf"

    plt.savefig(output_plot)

    print(f"[funcat] Dotplot saved to {output_plot}")

    return output_plot