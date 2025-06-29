import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: python bin_to_header.py <input_binary> <output_header> <variable_name>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    variable_name = sys.argv[3]

    try:
        with open(input_file, 'rb') as f_in:
            data = f_in.read()

        with open(output_file, 'w') as f_out:
            f_out.write(f'#pragma once\n\n')
            f_out.write(f'#include <cstddef>\n\n')
            f_out.write(f'const unsigned char {variable_name}[] = {{\n    ')
            
            for i, byte in enumerate(data):
                f_out.write(f'0x{byte:02x}, ')
                if (i + 1) % 16 == 0:
                    f_out.write('\n    ')
            
            f_out.write('\n};\n\n')
            f_out.write(f'const size_t {variable_name}_len = sizeof({variable_name});\n')

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
