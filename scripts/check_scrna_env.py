from importlib.metadata import version

import anndata as ad
import numpy as np
import scanpy as sc
import scvi
import torch


def main() -> None:
    print(f"scanpy: {version('scanpy')}")
    print(f"scvi-tools: {version('scvi-tools')}")
    print(f"torch: {torch.__version__}")
    print(f"cuda available: {torch.cuda.is_available()}")
    adata = ad.AnnData(np.ones((3, 4)))
    print(f"anndata test shape: {adata.shape}")


if __name__ == "__main__":
    main()
