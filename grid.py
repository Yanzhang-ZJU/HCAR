"""
MIT License

Copyright (c) 2023 Eternal (Shizhuo Cheng) @ Zhejiang University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

from biopandas.pdb import PandasPdb
from typing import Tuple, List
import numpy as np
import sys
import argparse
import logging
from pathlib import Path
from tqdm import tqdm


def load_pdb_to_dataframe(pdb_file: str) -> PandasPdb:
    """
    Load all atoms in a pdb file into a dataframe using biopandas.

    :param pdb_file: The path to the pdb file.
    :return: A dataframe containing all atoms in the pdb file.
    """
    ppdb = PandasPdb().read_pdb(pdb_file)
    return ppdb


def get_min_max_coordinates(
    ppdb: PandasPdb,
) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    """
    Get the min and max coordinates of all atoms on the x, y, and z axes.

    :param ppdb: The PandasPdb object containing all atoms.
    :return: A tuple containing the min and max coordinates on the x, y, and z axes.
    """
    df = ppdb.df["ATOM"]
    x_min, x_max = df["x_coord"].min(), df["x_coord"].max()
    y_min, y_max = df["y_coord"].min(), df["y_coord"].max()
    z_min, z_max = df["z_coord"].min(), df["z_coord"].max()

    return ((x_min, x_max), (y_min, y_max), (z_min, z_max))


def gen_grid(
    step: float,
    min_max_coordinates: Tuple[
        Tuple[float, float], Tuple[float, float], Tuple[float, float]
    ],
) -> List:
    """
    Generate a grid with each node as a list to store the index of df.

    :param step: The step size for the grid.
    :param min_max_coordinates: The min and max coordinates on the x, y, and z axes.
    :return: A grid with each node as a list.
    """
    x_range = np.arange(min_max_coordinates[0][0], min_max_coordinates[0][1], step)
    y_range = np.arange(min_max_coordinates[1][0], min_max_coordinates[1][1], step)
    z_range = np.arange(min_max_coordinates[2][0], min_max_coordinates[2][1], step)

    grid = [
        [[[] for _ in range(len(z_range))] for _ in range(len(y_range))]
        for _ in range(len(x_range))
    ]
    logging.info(f"xyz range: {min_max_coordinates}")
    logging.info(f"step size: {step}")
    
    return grid


def gen_bfgrid(grid: List[List[List[int]]]) -> np.ndarray:
    """
    Generate a b-factor grid and normalize it to [0, 100].

    :param grid: The grid with each node as a list to store the index of df.
    :return: A normalized b-factor grid.
    """
    # Generate a b-factor grid
    b_factor_grid = np.empty((len(grid), len(grid[0]), len(grid[0][0])))
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            for k in range(len(grid[0][0])):
                b_factor_grid[i][j][k] = len(grid[i][j][k])
                # Count the neighboring grids
                if i > 0:
                    b_factor_grid[i][j][k] += len(grid[i-1][j][k])
                if i < len(grid) - 1:
                    b_factor_grid[i][j][k] += len(grid[i+1][j][k])
                if j > 0:
                    b_factor_grid[i][j][k] += len(grid[i][j-1][k])
                if j < len(grid[0]) - 1:
                    b_factor_grid[i][j][k] += len(grid[i][j+1][k])
                if k > 0:
                    b_factor_grid[i][j][k] += len(grid[i][j][k-1])
                if k < len(grid[0][0]) - 1:
                    b_factor_grid[i][j][k] += len(grid[i][j][k+1])


    logging.info(f"Minimum b-factor: {b_factor_grid.min()}")
    logging.info(f"Maximum b-factor: {b_factor_grid.max()}")
    
    # Normalize the b-factor grid to [0, 100]
    b_factor_grid = (
        (b_factor_grid - b_factor_grid.min())
        / (b_factor_grid.max() - b_factor_grid.min())
        * 100
    )
    
    return b_factor_grid


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument("pdb_file", type=Path, help="The path to the pdb file.")
    parser.add_argument("step", type=float, help="The step size for the grid.")

    args = parser.parse_args()

    ppdb = load_pdb_to_dataframe(str(args.pdb_file))
    df = ppdb.df["ATOM"]
    bound = get_min_max_coordinates(ppdb)
    grid = gen_grid(args.step, bound)
    

    for index, atom in tqdm(df.iterrows()):
        x, y, z = atom["x_coord"], atom["y_coord"], atom["z_coord"]
        i = int((x - bound[0][0]) // args.step)
        j = int((y - bound[1][0]) // args.step)
        k = int((z - bound[2][0]) // args.step)
        grid[i][j][k].append(index)

    b_factor_grid = gen_bfgrid(grid)

    # Set the b-factor of atoms in the same grid node to this value in the df
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            for k in range(len(grid[0][0])):
                for index in grid[i][j][k]:
                    df.at[index, "b_factor"] = b_factor_grid[i][j][k]

    prefix = "new_"
    new_pdb_file = prefix + args.pdb_file.name
    ppdb.to_pdb(path=str(new_pdb_file), records=["ATOM"])


if __name__ == "__main__":
    main()

