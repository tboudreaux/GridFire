# generate_partition_binary.py
import sys
import struct
import re

def generate_binary_file(input_filepath, output_filepath):
    RECORD_FORMAT = '<i i d 24d'

    record_regex = re.compile(r"""
        ^\s*[a-zA-Z0-9]+\s*\n
        ^\s*
        (?P<z>\d+)\s+
        (?P<a>\d+)\s+
        (?P<spin>[\d\.]+)
        \s*\n
        (?P<coeffs>
          (?:^\s*(?:[\d.]+\s+){7}[\d.]+)\s*\n
          (?:^\s*(?:[\d.]+\s+){7}[\d.]+)\s*\n
          (?:^\s*(?:[\d.]+\s+){7}[\d.]+)
        )
    """, re.MULTILINE | re.VERBOSE)

    with open(input_filepath, 'r') as f:
        content = f.read()

    record_count = 0
    with open(output_filepath, 'wb') as f_out:
        for match in record_regex.finditer(content):
            z = int(match.group('z'))
            a = int(match.group('a'))
            spin = float(match.group('spin'))
            coeffs_str = match.group('coeffs')
            g_values = [float(val) for val in coeffs_str.split()]
            print(f"Processing Z={z}, A={a}, Spin={spin}")

            if len(g_values) != 24:
                print(f"Warning: Found {len(g_values)} coefficients for Z={z}, A={a}. Expected 24. Skipping.")
                continue

            packed_data = struct.pack(RECORD_FORMAT, z, a, spin, *g_values)
            f_out.write(packed_data)
            record_count += 1

    print(f"Found and processed {record_count} records. Wrote to binary file {output_filepath}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input_partition_file> <output_binary_file>")
        sys.exit(1)
    generate_binary_file(sys.argv[1], sys.argv[2])