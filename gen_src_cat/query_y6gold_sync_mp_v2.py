import pyvo as vo
import time
import os
import numpy as np
from multiprocessing import Pool

# --- Configuration ---
OUT_DIR = "des_y6_chunks"
CONCURRENT_PROCESSES = 4  # Reduce to 2-3 if the TAP service throttles requests.
TAP_URL = "https://datalab.noirlab.edu/tap"
MAX_ROWS = 300000         # Maximum number of rows requested by each synchronous query.


# ==========================================
# Function: format_dec
# Method: Convert declination to the pXX/mXX convention used by the pipeline.
# ==========================================
def format_dec(dec_str):
    """
    Format a declination boundary for a catalog tile name.

    Discard the fractional part, prefix non-negative values with ``p``,
    replace the minus sign with ``m`` for negative values, and zero-pad the
    absolute value to two digits (for example, ``p05`` or ``m09``).
    """
    dec_int = int(float(dec_str))

    if dec_int >= 0:
        return f"p{dec_int:02d}"
    return f"m{abs(dec_int):02d}"


# ==========================================
# Function: build_output_filename
# Method: Build the tile name expected by the Fourier_Quad pipeline.
# ==========================================
def build_output_filename(ra1, ra2, dec1, dec2):
    """Return a tile name with padded RA and signed Dec boundaries."""
    ra1_fmt = f"{int(float(ra1)):03d}"
    ra2_fmt = f"{int(float(ra2)):03d}"
    dec1_fmt = format_dec(dec1)
    dec2_fmt = format_dec(dec2)
    return f"des_y6_RA_{ra1_fmt}_{ra2_fmt}_Dec_{dec1_fmt}_{dec2_fmt}.dat"


# ==========================================
# Function: get_file_row_count
# Method: Count data rows in an existing ASCII table without loading it.
# ==========================================
def get_file_row_count(filepath):
    """Return the number of data rows after subtracting the header line."""
    try:
        with open(filepath, "r") as f:
            # Streaming the file avoids loading a potentially large tile into memory.
            return sum(1 for _ in f) - 1  # Subtract the header row.
    except:
        return 0


# ==========================================
# Function: download_one_chunk_sync
# Method: Run one synchronous TAP/ADQL sky-tile query per worker process.
# ==========================================
def download_one_chunk_sync(coords):
    """
    Download one rectangular sky tile with a synchronous query.

    ``coords`` contains ``(ra1, ra2, dec1, dec2)`` in degrees.
    """
    ra1, ra2, dec1, dec2 = coords
    filename = build_output_filename(ra1, ra2, dec1, dec2)
    output_file = os.path.join(OUT_DIR, filename)

    if os.path.exists(output_file):
        count = get_file_row_count(output_file)
        if count >= MAX_ROWS:
            print(
                f"[!] Row-limit tile detected ({count} rows): {filename}. "
                "Deleting it before retrying the same tile...",
                flush=True,
            )
            os.remove(output_file)
            # The current implementation retries the same tile; it does not
            # automatically subdivide a truncated query.
        else:
            # print(f"[-] Skipping complete tile ({count} rows): {filename}", flush=True)
            return

    # Optional simple resume check (superseded by the row-count check above).
    # if os.path.exists(output_file):
    #     print(f"[-] Skipping existing tile: {filename}", flush=True)
    #     return

    # Create a TAP service connection dedicated to this worker process.
    try:
        service = vo.dal.TAPService(TAP_URL)

        # Prepare cone geometry for an optional spatial-query implementation.
        # The rectangular ADQL query below currently uses only the tile bounds.
        ra_center = (ra1 + ra2) / 2.0
        dec_center = (dec1 + dec2) / 2.0
        radius = 0.75  # A 1 x 1 degree tile fits inside a 0.75 degree cone.

        # Build a rectangular ADQL query with an explicit row limit.
        adql_query = f"""
            SELECT TOP {MAX_ROWS}
                flags_footprint, flags_foreground, flags_gold, ext_mash,
                ra, dec, bdf_mag_g, bdf_mag_err_g, bdf_mag_r, bdf_mag_err_r,
                bdf_mag_i, bdf_mag_err_i, bdf_mag_z, bdf_mag_err_z, bdf_mag_y, bdf_mag_err_y,
                dnf_z, dnf_zsigma
            FROM
                des_dr2.y6_gold
            WHERE
                ra >= {ra1} AND ra < {ra2}
                AND dec >= {dec1} AND dec < {dec2}
            """

        # print(f"[+] Starting synchronous query: {filename}...", flush=True)

        # Run the synchronous query. TAPService.search blocks until the result
        # is downloaded or the service reports an error or timeout.
        start_time = time.time()
        result = service.search(adql_query)
        elapsed_time = time.time() - start_time

        # Convert the result and write the pipeline-compatible ASCII table.
        results_table = result.to_table()
        num_rows = len(results_table)

        if num_rows == 0:
            # print(f"[!] Empty tile skipped: {filename} ({elapsed_time:.1f}s)", flush=True)
            return

        results_table.write(output_file, format="ascii.commented_header", overwrite=True)

        # Warn when TOP may have truncated the result.
        if num_rows >= MAX_ROWS:
            print(
                f"[!!!] WARNING: {filename} reached MAX_ROWS ({MAX_ROWS}); "
                "the tile may be incomplete.",
                flush=True,
            )
        else:
            print(
                f"[OK] Downloaded: {filename} ({num_rows} rows, "
                f"elapsed: {elapsed_time:.1f}s)",
                flush=True,
            )

        time.sleep(3)

    except Exception as e:
        # Synchronous queries can encounter 504 Gateway Timeout or ReadTimeout errors.
        print(
            f"[ERROR] Worker failed while processing {filename} "
            f"(possible network timeout or server rejection): {e}",
            flush=True,
        )

        time.sleep(10)


# ==========================================
# Function: main
# Method: Build a 1-degree grid and download its tiles with a process pool.
# ==========================================
def main():
    """Create the output directory and dispatch all configured sky tiles."""

    # Create the output directory.
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    # Configure the target sky range in degrees.
    TARGET_RA_MIN, TARGET_RA_MAX = 299, 360
    TARGET_DEC_MIN, TARGET_DEC_MAX = -80.0, 20.0
    CHUNK_SIZE = 1.0

    # Build the list of rectangular tile queries.
    ra_bins = np.arange(TARGET_RA_MIN, TARGET_RA_MAX, CHUNK_SIZE)
    dec_bins = np.arange(TARGET_DEC_MIN, TARGET_DEC_MAX, CHUNK_SIZE)

    task_list = []
    for ra1 in ra_bins:
        for dec1 in dec_bins:
            task_list.append((ra1, ra1 + CHUNK_SIZE, dec1, dec1 + CHUNK_SIZE))

    print(f"Total sky-tile tasks: {len(task_list)}")
    print(f"Worker processes: {CONCURRENT_PROCESSES} (synchronous query mode)")

    # Dispatch the synchronous queries through the process pool.
    with Pool(processes=CONCURRENT_PROCESSES) as pool:
        pool.map(download_one_chunk_sync, task_list)


if __name__ == "__main__":
    main()
