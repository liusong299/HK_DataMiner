import typer
from hkdataminer import workflows
import os
import subprocess
from pathlib import Path

app = typer.Typer(help="HKDataMiner CLI")
cluster_app = typer.Typer(help="Clustering commands")
lump_app = typer.Typer(help="Lumping commands")
app.add_typer(cluster_app, name="cluster")
app.add_typer(lump_app, name="lump")

@cluster_app.command("kcenters")
def kcenters(
    trajlist: str = typer.Option("trajlist", help="List of trajectory files"),
    atomlist: str = typer.Option("atom_indices", help="List of atom index files"),
    topology: str = typer.Option("native.pdb", help="Topology file"),
    homedir: str = typer.Option(".", help="Home directory"),
    iext: str = typer.Option("xtc", help="Input extension"),
    n_clusters: int = typer.Option(100, help="Number of clusters"),
    stride: int = typer.Option(None, help="Stride"),
    output_dir: str = typer.Option(".", help="Output directory"),
):
    """Run K-Centers clustering."""
    workflows.run_clustering(
        trajListFns=trajlist,
        atomListFns=atomlist,
        topology=topology,
        homedir=homedir,
        iext=iext,
        n_clusters=n_clusters,
        stride=stride,
        output_dir=output_dir
    )

@cluster_app.command("aplod")
def aplod(
    trajlist: str = typer.Option("trajlist", help="List of trajectory files"),
    atomlist: str = typer.Option("atom_indices", help="List of atom index files"),
    topology: str = typer.Option("native.pdb", help="Topology file"),
    homedir: str = typer.Option(".", help="Home directory"),
    iext: str = typer.Option("xtc", help="Input extension"),
    rho_cutoff: float = typer.Option(1.0, help="Rho cutoff"),
    delta_cutoff: float = typer.Option(1.0, help="Delta cutoff"),
    n_neighbors: int = typer.Option(100, help="Number of neighbors"),
    stride: int = typer.Option(None, help="Stride"),
    output_dir: str = typer.Option(".", help="Output directory"),
):
    """Run APLoD clustering (on Phi/Psi angles)."""
    workflows.run_aplod(
        trajListFns=trajlist,
        atomListFns=atomlist,
        topology=topology,
        homedir=homedir,
        iext=iext,
        rho_cutoff=rho_cutoff,
        delta_cutoff=delta_cutoff,
        n_neighbors=n_neighbors,
        stride=stride,
        output_dir=output_dir
    )

@lump_app.command("pcca")
def pcca(
    assignments: str = typer.Option(..., help="Assignments file"),
    traj_len: str = typer.Option("traj_len.txt", help="Trajectory lengths file"),
    n_macro: int = typer.Option(6, help="Number of macrostates"),
    homedir: str = typer.Option(".", help="Output directory"),
):
    """Run PCCA lumping."""
    workflows.run_lumping(
        assignments_file=assignments,
        traj_len_file=traj_len,
        n_macro_states=n_macro,
        homedir=homedir
    )

@app.command("tutorial")
def tutorial():
    """Run the end-to-end tutorial."""
    script_path = Path("scripts/run_tutorial.sh")
    if not script_path.exists():
        # Fallback if running from installed package and scripts not present?
        # But we are in the repo.
        typer.echo("scripts/run_tutorial.sh not found. Are you in the project root?")
        raise typer.Exit(1)
    
    typer.echo("Running tutorial script...")
    subprocess.run(["bash", str(script_path)], check=True)

if __name__ == "__main__":
    app()