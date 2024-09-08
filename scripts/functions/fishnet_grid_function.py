import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import shapely
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

def process_polygon(row, spatial_resolution, crs):
    buffer_distance = 0.001  # adjust this value as needed

    xmin, ymin, xmax, ymax = row['geometry'].bounds
    xmin -= buffer_distance
    ymin -= buffer_distance
    xmax += buffer_distance
    ymax += buffer_distance
    
    cell_size = spatial_resolution
    grid_cells = []

    for x0 in np.arange(xmin, xmax, cell_size):
        for y0 in np.arange(ymin, ymax, cell_size):
            x1 = x0 + cell_size
            y1 = y0 + cell_size
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))

    return grid_cells

def create_grids_parallel(
    df,
    longitude_column,
    latitude_column,
    csv_crs,
    grid_crs,
    transformation_crs,
    final_cell_size,
    coarse_cell_size=None,
    parallel=True
):
    # Convert latitude and longitude columns into Point geometries
    geometry = [Point(xy) for xy in zip(df[longitude_column], df[latitude_column])]

    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs=csv_crs)

    # First creating a coarse grid
    buffer_distance = 0.001  # adjust this value as needed

    xmin, ymin, xmax, ymax = gdf.total_bounds
    xmin -= buffer_distance
    ymin -= buffer_distance
    xmax += buffer_distance
    ymax += buffer_distance


    # Determine cell size and number of cells
    n_cells = 50
    if coarse_cell_size is None:
        coarse_cell_size = (xmax - xmin) / n_cells
    else:
        n_cells = int((xmax - xmin) / coarse_cell_size)

    # Create coarse grid cells
    grid_cells = []
    for x0 in np.arange(xmin, xmax + coarse_cell_size, coarse_cell_size):
        for y0 in np.arange(ymin, ymax + coarse_cell_size, coarse_cell_size):
            x1 = x0 + coarse_cell_size
            y1 = y0 + coarse_cell_size
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))

    coarse_grid = gpd.GeoDataFrame(grid_cells, columns=['geometry'], crs=grid_crs)

    # Filter out unused grids
    filtered_polygons = gpd.sjoin(coarse_grid, gdf, how="inner", predicate='intersects').drop_duplicates(subset='geometry')

    spatial_resolution = final_cell_size

    if parallel:
        # Parallel processing to create fine grid
        grid_cells = Parallel(n_jobs=-1)(
            delayed(process_polygon)(row, spatial_resolution, filtered_polygons.crs)
            for _, row in tqdm(filtered_polygons.iterrows(), total=len(filtered_polygons), desc="Processing Polygons", position=0, leave=True)
        )
    else:
        # Sequential processing to create fine grid
        grid_cells = [
            process_polygon(row, spatial_resolution, filtered_polygons.crs)
            for _, row in tqdm(filtered_polygons.iterrows(), total=len(filtered_polygons), desc="Processing Polygons", position=0, leave=True)
        ]

    # Flatten list of grid cells
    grid_cells_flat = [item for sublist in grid_cells for item in sublist]

    fine_grid = gpd.GeoDataFrame(grid_cells_flat, columns=['geometry'], crs=transformation_crs)

    # Filter out unused grids again
    filtered_polygons_fine = gpd.sjoin(fine_grid, gdf, how="inner", predicate='intersects').drop_duplicates(subset='geometry')

    filtered_polygons_fine = filtered_polygons_fine['geometry']
    return filtered_polygons_fine
