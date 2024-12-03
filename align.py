import re
import argparse

def align_indentation(file_path, output_path, indent_width=4):
    """
    Align indentation for C++ code in the given file, ignoring lines starting with "//".
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    indent_level = 0
    indent_space = ' ' * indent_width
    formatted_lines = []

    for line in lines:
        stripped = line.strip()

        # Ignore single-line comments
        if stripped.startswith('//'):
            formatted_lines.append(line)
            continue

        # Skip empty lines
        if not stripped:
            formatted_lines.append(line)
            continue

        # Adjust indent level based on braces
        if stripped.endswith('}') or stripped.startswith('}'):
            indent_level -= 1

        # Apply the current indent
        formatted_line = f"{indent_space * indent_level}{stripped}\n"
        formatted_lines.append(formatted_line)

        # Adjust indent level for opening braces
        if stripped.endswith('{') or '{' in stripped and not stripped.endswith(';'):
            indent_level += 1

        # Prevent negative indentation
        indent_level = max(indent_level, 0)

    # Write to the output file
    with open(output_path, 'w') as output_file:
        output_file.writelines(formatted_lines)

    print(f"Formatted code written to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align indentation for C++ code.")
    parser.add_argument("input_file", help="Path to the input C++ file.")
    parser.add_argument("output_file", help="Path to the output formatted file.")
    parser.add_argument("--indent-width", type=int, default=4, help="Number of spaces for each indentation level.")

    args = parser.parse_args()

    align_indentation(args.input_file, args.output_file, args.indent_width)



    
