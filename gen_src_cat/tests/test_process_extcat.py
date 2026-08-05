#!/usr/bin/env python3

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


CANONICAL_COLUMNS = [
    "flags_footprint",
    "flags_foreground",
    "flags_gold",
    "ext_mash",
    "ra",
    "dec",
    "bdf_mag_g",
    "bdf_mag_err_g",
    "bdf_mag_r",
    "bdf_mag_err_r",
    "bdf_mag_i",
    "bdf_mag_err_i",
    "bdf_mag_z",
    "bdf_mag_err_z",
    "bdf_mag_y",
    "bdf_mag_err_y",
    "dnf_z",
    "dnf_zsigma",
]


# ==========================================
# Function: Build one canonical synthetic catalog row
# Method: Use distinct finite numeric strings while injecting selected celestial coordinates.
# ==========================================
def canonical_row(seed, ra, dec):
    return [
        str(seed),
        str(seed + 1),
        str(seed + 2),
        str(seed + 3),
        str(ra),
        str(dec),
        f"{20.0 + seed / 10.0:.3f}",
        "0.010",
        f"{19.0 + seed / 10.0:.3f}",
        "0.011",
        f"{18.0 + seed / 10.0:.3f}",
        "0.012",
        f"{17.0 + seed / 10.0:.3f}",
        "0.013",
        f"{16.0 + seed / 10.0:.3f}",
        "0.014",
        f"{0.5 + seed / 100.0:.3f}",
        "0.050",
    ]


# ==========================================
# Function: Return the exact canonical output header
# Method: Join the schema names with the Astropy-compatible leading hash marker.
# ==========================================
def canonical_header():
    return "# " + " ".join(CANONICAL_COLUMNS)


# ==========================================
# Function: Format a declination tile boundary
# Method: Match the pXX/mXX convention used by the C++ and Python generators.
# ==========================================
def format_dec(value):
    return f"{'p' if value >= 0 else 'm'}{abs(value):02d}"


# ==========================================
# Function: Build one expected output tile basename
# Method: Apply three-digit RA padding and signed declination boundaries.
# ==========================================
def tile_name(ra_lower, dec_lower):
    return (
        f"des_y6_RA_{ra_lower:03d}_{ra_lower + 1:03d}_Dec_"
        f"{format_dec(dec_lower)}_{format_dec(dec_lower + 1)}.dat"
    )


# ==========================================
# Function: Run the standalone MPI test executable
# Method: Use an oversubscribed local communicator and preserve diagnostics on failure.
# ==========================================
def run_process(executable, ranks, arguments, expect_success=True):
    mpirun = os.environ.get("MPIRUN", shutil.which("mpirun") or "mpirun")
    command = [mpirun, "--oversubscribe", "-np", str(ranks), str(executable), *arguments]
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    if expect_success and completed.returncode != 0:
        raise AssertionError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    if not expect_success and completed.returncode == 0:
        raise AssertionError(f"command unexpectedly succeeded: {' '.join(command)}")
    return completed


# ==========================================
# Function: Read every generated tile as an ordered byte mapping
# Method: Select only canonical tile basenames and preserve exact contents for comparison.
# ==========================================
def read_tiles(directory):
    return {
        path.name: path.read_bytes()
        for path in sorted(directory.glob("des_y6_RA_*_Dec_*.dat"))
    }


# ==========================================
# Function: Write mixed-delimiter inputs with one shared schema
# Method: Keep identical column order across CSV and whitespace files, add one malformed row,
#         and place one selected file in a nested directory.
# ==========================================
def write_mixed_inputs(input_directory):
    input_directory.mkdir(parents=True)
    nested = input_directory / "nested"
    nested.mkdir()

    csv_order = list(range(len(CANONICAL_COLUMNS)))
    row_a = canonical_row(10, 299.2, -79.2)
    row_b = canonical_row(20, 359.9, 5.4)
    with (input_directory / "keep_a.csv").open("w", encoding="utf-8", newline="\n") as output:
        output.write("\ufeff" + ",".join(CANONICAL_COLUMNS[index] for index in csv_order) + "\n")
        output.write(",".join(row_a[index] for index in csv_order) + "\n")
        output.write("this,is,a,malformed,row\n")
        output.write(",".join(row_b[index] for index in csv_order) + "\n")

    row_c = canonical_row(30, 360.0, 90.0)
    row_d = canonical_row(40, 0.1, -0.1)
    with (nested / "keep_b.cat").open("w", encoding="utf-8", newline="\n") as output:
        output.write(canonical_header() + "\n")
        output.write(" ".join(row_c) + "\n")
        output.write(" ".join(row_d) + "\n")

    ignored = canonical_row(90, 100.0, 10.0)
    with (input_directory / "ignore.cat").open("w", encoding="utf-8", newline="\n") as output:
        output.write(canonical_header() + "\n")
        output.write(" ".join(ignored) + "\n")

    return {
        tile_name(299, -80): row_a,
        tile_name(359, 5): row_b,
        tile_name(0, 89): row_c,
        tile_name(0, -1): row_d,
    }


# ==========================================
# Function: Assert exact mixed-input tile contents
# Method: Verify one canonical header and one projected row in each expected sky tile.
# ==========================================
def assert_mixed_outputs(output_directory, expected_rows):
    tiles = read_tiles(output_directory)
    if set(tiles) != set(expected_rows):
        raise AssertionError(f"unexpected tile set: {sorted(tiles)}")
    for name, row in expected_rows.items():
        expected = (canonical_header() + "\n" + " ".join(row) + "\n").encode()
        if tiles[name] != expected:
            raise AssertionError(f"unexpected contents for {name}: {tiles[name]!r}")


# ==========================================
# Function: Test pass-through columns and MPI determinism
# Method: Run identical mixed-delimiter schemas with one and three ranks, compare bytes, and
#         verify filename filtering, recursion, malformed skipping, and pole boundaries.
# ==========================================
def test_parallel_projection(executable, workspace):
    input_directory = workspace / "mixed_input"
    expected_rows = write_mixed_inputs(input_directory)
    output_one = workspace / "mixed_output_one"
    output_three = workspace / "mixed_output_three"
    arguments = [
        "--input-dir",
        str(input_directory),
        "--contains",
        "keep",
        "--recursive",
        "true",
        "--delimiter",
        "auto",
        "--header",
        "auto",
        "--malformed",
        "skip",
    ]
    failed_output = workspace / "mixed_output_fail"
    run_process(
        executable,
        2,
        [*arguments, "--malformed", "fail", "--output-dir", str(failed_output)],
        expect_success=False,
    )
    if failed_output.exists() and read_tiles(failed_output):
        raise AssertionError("fail-on-malformed published final tiles")
    if failed_output.exists() and list(failed_output.glob(".process_extcat_staging_*")):
        raise AssertionError("fail-on-malformed left temporary staging data")

    run_process(executable, 1, [*arguments, "--output-dir", str(output_one)])
    run_process(executable, 3, [*arguments, "--output-dir", str(output_three)])
    assert_mixed_outputs(output_one, expected_rows)
    assert_mixed_outputs(output_three, expected_rows)
    if read_tiles(output_one) != read_tiles(output_three):
        raise AssertionError("one-rank and three-rank outputs differ")


# ==========================================
# Function: Test arbitrary ordered projection and output lifecycle policies
# Method: Apply the requested {5,3,4,1} selection to a seven-field table, retain a string
#         payload, then verify nested-output rejection and transactional overwrite behavior.
# ==========================================
def test_explicit_columns_and_overwrite(executable, workspace):
    input_directory = workspace / "explicit_input"
    input_directory.mkdir()
    header = ["c1", "c2", "c3", "c4", "ra", "dec", "tail"]
    row = ["source-1", "20", "30", "40", "42.25", "-10.75", "tail-value"]
    raw_path = input_directory / "projection_raw.txt"
    raw_path.write_text(
        "# " + " ".join(header) + "\n" + " ".join(row) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    output_directory = workspace / "explicit_output"
    columns = "5,3,4,1"
    arguments = [
        "--input-dir",
        str(input_directory),
        "--output-dir",
        str(output_directory),
        "--contains",
        "projection",
        "--header",
        "present",
        "--columns",
        columns,
    ]
    run_process(
        executable,
        2,
        [
            "--input-dir",
            str(input_directory),
            "--output-dir",
            str(input_directory / "unsafe_nested_output"),
            "--contains",
            "projection",
            "--header",
            "present",
            "--columns",
            columns,
        ],
        expect_success=False,
    )
    run_process(executable, 2, arguments)

    expected_name = tile_name(42, -11)
    expected_bytes = ("# ra c3 c4 c1\n42.25 30 40 source-1\n").encode()
    if read_tiles(output_directory) != {expected_name: expected_bytes}:
        raise AssertionError("explicit column projection produced unexpected output")

    run_process(executable, 2, arguments, expect_success=False)
    stale = output_directory / tile_name(123, 0)
    stale.write_text(canonical_header() + "\n", encoding="utf-8", newline="\n")
    note = output_directory / "preserve_me.txt"
    note.write_text("unrelated\n", encoding="utf-8", newline="\n")
    run_process(executable, 2, [*arguments, "--existing", "overwrite"])
    if stale.exists():
        raise AssertionError("overwrite did not remove a stale generated tile")
    if note.read_text(encoding="utf-8") != "unrelated\n":
        raise AssertionError("overwrite changed an unrelated output file")
    if read_tiles(output_directory) != {expected_name: expected_bytes}:
        raise AssertionError("overwrite did not publish the expected tile set")


# ==========================================
# Function: Test variable-width pass-through output
# Method: Leave explicit projection disabled and verify that all six named input fields,
#         including string payloads, retain their original order.
# ==========================================
def test_variable_width_pass_through(executable, workspace):
    input_directory = workspace / "passthrough_input"
    input_directory.mkdir()
    raw_path = input_directory / "passthrough.cat"
    header = "# object_id ra dec flux note extra"
    row = "obj-A 81.25 1.75 99.5 keep-me final"
    raw_path.write_text(header + "\n" + row + "\n", encoding="utf-8", newline="\n")

    output_directory = workspace / "passthrough_output"
    run_process(
        executable,
        2,
        [
            "--input-dir",
            str(input_directory),
            "--output-dir",
            str(output_directory),
            "--contains",
            "passthrough",
        ],
    )
    expected = (header + "\n" + row + "\n").encode()
    if read_tiles(output_directory) != {tile_name(81, 1): expected}:
        raise AssertionError("disabled projection did not preserve all input columns")


# ==========================================
# Function: Test headerless coordinate/output decoupling
# Method: Supply raw RA/Dec indices independently from a two-column output projection and
#         verify deterministic generic header names.
# ==========================================
def test_headerless_coordinate_columns(executable, workspace):
    input_directory = workspace / "headerless_input"
    input_directory.mkdir()
    raw_path = input_directory / "headerless.cat"
    raw_path.write_text(
        "source-A 12.25 keep -1.25 9\n",
        encoding="utf-8",
        newline="\n",
    )

    output_directory = workspace / "headerless_output"
    run_process(
        executable,
        3,
        [
            "--input-dir",
            str(input_directory),
            "--output-dir",
            str(output_directory),
            "--header",
            "absent",
            "--columns",
            "5,1",
            "--ra-column",
            "2",
            "--dec-column",
            "4",
        ],
    )
    expected = b"# column_5 column_1\n9 source-A\n"
    if read_tiles(output_directory) != {tile_name(12, -2): expected}:
        raise AssertionError("headerless coordinate columns were coupled to output order")


# ==========================================
# Function: Test incompatible pass-through schemas fail before publication
# Method: Present two different output headers and ensure no mixed-schema tile is committed.
# ==========================================
def test_schema_mismatch_rejected(executable, workspace):
    input_directory = workspace / "schema_mismatch_input"
    input_directory.mkdir()
    (input_directory / "schema_a.cat").write_text(
        "# id ra dec\na 20.0 2.0\n", encoding="utf-8", newline="\n"
    )
    (input_directory / "schema_b.cat").write_text(
        "# id ra dec flux\nb 20.1 2.1 5\n", encoding="utf-8", newline="\n"
    )
    output_directory = workspace / "schema_mismatch_output"
    run_process(
        executable,
        2,
        [
            "--input-dir",
            str(input_directory),
            "--output-dir",
            str(output_directory),
            "--contains",
            "schema_",
        ],
        expect_success=False,
    )
    if output_directory.exists() and read_tiles(output_directory):
        raise AssertionError("incompatible schemas published output tiles")


# ==========================================
# Function: Run all dependency-free process_extcat integration tests
# Method: Resolve the executable, isolate fixtures in a temporary directory, and fail fast.
# ==========================================
def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_process_extcat.py /path/to/process_extcat")
    executable = Path(sys.argv[1]).resolve()
    if not executable.is_file():
        raise SystemExit(f"test executable does not exist: {executable}")

    with tempfile.TemporaryDirectory(prefix="process_extcat_test_") as temporary:
        workspace = Path(temporary)
        test_parallel_projection(executable, workspace)
        test_explicit_columns_and_overwrite(executable, workspace)
        test_variable_width_pass_through(executable, workspace)
        test_headerless_coordinate_columns(executable, workspace)
        test_schema_mismatch_rejected(executable, workspace)
    print("process_extcat integration tests: PASS")


if __name__ == "__main__":
    main()
