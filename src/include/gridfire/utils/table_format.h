#pragma once

#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <memory>
#include <print>
#include <cwchar>

#include "nlohmann/json.hpp"



namespace gridfire::utils {
    inline size_t visual_width(const std::string& s) {
        // IMPORTANT: std::setlocale(LC_ALL, "") must be called once in main()
        // for mbtowc and wcwidth to function correctly with the system's locale.

        size_t width = 0;
        std::mbtowc(nullptr, nullptr, 0); // Reset multi-byte state

        const char* p = s.c_str();
        const char* end = s.c_str() + s.length();

        while (p < end) {
            wchar_t wc;
            // Convert the next multibyte char to a wide char
            const int byte_len = std::mbtowc(&wc, p, end - p);

            if (byte_len <= 0) {
                // Invalid byte sequence or null char.
                // Treat as a 1-width '?' and advance by 1 byte to avoid infinite loop.
                width++;
                p++;
                continue;
            }

            // Get the visual width of the wide char
            int char_width = wcwidth(wc);
            if (char_width != -1) {
                width += char_width;
            }
            // else: char_width == -1 means non-printable/control char; treat as 0 width

            p += byte_len; // Advance by the number of bytes consumed
        }
        return width;
    }

    class ColumnBase {
    public:
        virtual ~ColumnBase() = default;
        // Gets the string representation of the data at a given row
        [[nodiscard]] virtual std::string getCellData(size_t rowIndex) const = 0;
        // Gets the header text for the column
        [[nodiscard]] virtual std::string getHeader() const = 0;
        // Gets the number of data rows in the column
        [[nodiscard]] virtual size_t getRowCount() const = 0;
    };

    template<typename T>
    class Column final : public ColumnBase {
    public:
        Column(std::string header, const std::vector<T>& data)
            : m_header(std::move(header)), m_data(data) {}

        [[nodiscard]] std::string getCellData(size_t rowIndex) const override {
            std::stringstream ss;
            if (rowIndex < m_data.size()) {
                ss << m_data[rowIndex];
            }
            return ss.str();
        }

        [[nodiscard]] std::string getHeader() const override {
            return m_header;
        }

        [[nodiscard]] size_t getRowCount() const override {
            return m_data.size();
        }
    private:
        std::string m_header;
        std::vector<T> m_data;

    };

    class TableBase {
    public:
        virtual ~TableBase() = default;
        virtual size_t ncols() const = 0;
        virtual size_t nrows() const = 0;
        virtual const ColumnBase& operator[](size_t i) const = 0;
        virtual const std::vector<std::unique_ptr<ColumnBase>> getColumns() const;
        virtual std::string getName() const = 0;
        virtual void toCSV(const std::string& filename) const = 0;
    };

    template <typename T>
    class Table final : public TableBase {
    public:
        Table(const std::string& name, const std::vector<Column<T>>& columns)
            : m_name(name) {
            for (const auto& col : columns) {
                m_columns.push_back(std::make_unique<Column<T>>(col));
            }
        }

        Table(const std::string& name, const std::unordered_map<std::string, std::vector<T>>& columnData)
            : m_name(name) {
            for (const auto& [header, data] : columnData) {
                m_columns.push_back(std::make_unique<Column<T>>(header, data));
            }
        }

        Table(const std::string& name, const std::map<std::string, std::vector<T>>& columnData)
            : m_name(name) {
            for (const auto& [header, data] : columnData) {
                m_columns.push_back(std::make_unique<Column<T>>(header, data));
            }
        }

        Table(const std::string& name, const std::unordered_map<fourdst::atomic::Species, std::vector<T>>& columnData)
            : m_name(name) {
            for (const auto& [species, data] : columnData) {
                m_columns.push_back(std::make_unique<Column<T>>(std::string(species.name()), data));
            }
        }

        Table(const std::string& name, const std::map<fourdst::atomic::Species, std::vector<T>>& columnData)
            : m_name(name) {
            for (const auto& [species, data] : columnData) {
                m_columns.push_back(std::make_unique<Column<T>>(std::string(species.name()), data));
            }
        }

        size_t ncols() const override {
            return m_columns.size();
        }

        size_t nrows() const override {
            return m_columns.empty() ? 0 : m_columns[0]->getRowCount();
        }

        const ColumnBase& operator[](size_t i) const override {
            return *(m_columns[i]);
        }

        const std::vector<std::unique_ptr<ColumnBase>> getColumns() const override {
            std::vector<std::unique_ptr<ColumnBase>> cols;
            for (const auto& col : m_columns) {
                cols.push_back(std::make_unique<Column<T>>(*col));
            }
            return cols;
        }

        std::string getName() const override {
            return m_name;
        }

        void toCSV(const std::string& filename) const override {
            std::ofstream csvFile(filename);
            if (!csvFile.is_open()) {
                throw std::runtime_error("Failed to open file for writing: " + filename);
            }

            // Write header
            for (size_t j = 0; j < m_columns.size(); ++j) {
                csvFile << m_columns[j]->getHeader();
                if (j < m_columns.size() - 1) {
                    csvFile << ",";
                }
            }
            csvFile << "\n";

            // Write data rows
            size_t num_rows = nrows();
            for (size_t i = 0; i < num_rows; ++i) {
                for (size_t j = 0; j < m_columns.size(); ++j) {
                    csvFile << m_columns[j]->getCellData(i);
                    if (j < m_columns.size() - 1) {
                        csvFile << ",";
                    }
                }
                csvFile << "\n";
            }

            csvFile.close();
        }
    private:
        std::string m_name;
        std::vector<std::unique_ptr<Column<T>>> m_columns;
    };


    inline std::string format_table(const std::string& tableName, const std::vector<std::unique_ptr<ColumnBase>>& columns) {
        // --- 1. Handle Empty Table ---
        if (columns.empty()) {
            return tableName + "\n(Table has no columns)\n";
        }

        // --- 2. Determine dimensions and calculate column widths ---
        size_t num_cols = columns.size();
        size_t num_rows = 0;
        for(const auto& col : columns) {
            num_rows = std::max(num_rows, col->getRowCount());
        }

        std::vector<size_t> col_widths(num_cols);
        for (size_t j = 0; j < num_cols; ++j) {
            col_widths[j] = columns[j]->getHeader().length();
            for (size_t i = 0; i < num_rows; ++i) {
                col_widths[j] = std::max(col_widths[j], columns[j]->getCellData(i).length());
            }
        }

        // --- 3. Build the table string using stringstream ---
        std::stringstream table_ss;

        // --- Table Title ---
        const size_t total_width = std::accumulate(col_widths.begin(), col_widths.end(), 0) + (num_cols * 3) + 1; // NOLINT(*-fold-init-type)
        const size_t title_padding = (total_width > tableName.length()) ? (total_width - tableName.length()) / 2 : 0;
        table_ss << std::string(title_padding, ' ') << tableName << "\n";

        // --- Helper to draw horizontal border ---
        auto draw_border = [&]() {
            table_ss << "+";
            for (const size_t width : col_widths) {
                table_ss << std::string(width + 2, '-'); // +2 for padding
                table_ss << "+";
            }
            table_ss << "\n";
        };

        // --- Draw Top Border ---
        draw_border();

        // --- Draw Header Row ---
        table_ss << "|";
        for (size_t j = 0; j < num_cols; ++j) {
            table_ss << " " << std::left << std::setw(col_widths[j]) << columns[j]->getHeader() << " |"; // NOLINT(*-narrowing-conversions)
        }
        table_ss << "\n";

        // --- Draw Separator ---
        draw_border();

        // --- Draw Data Rows ---
        for (size_t i = 0; i < num_rows; ++i) {
            table_ss << "|";
            for (size_t j = 0; j < num_cols; ++j) {
                table_ss << " " << std::left << std::setw(col_widths[j]) << columns[j]->getCellData(i) << " |"; // NOLINT(*-narrowing-conversions)
            }
            table_ss << "\n";
        }

        // --- Draw Bottom Border ---
        draw_border();

        return table_ss.str();
    }

    inline std::string format_table(const TableBase& table) {
        return format_table(table.getName(), table.getColumns());
    }

    inline void print_table(const std::string& tableName, const std::vector<std::unique_ptr<ColumnBase>>& columns) {
        // --- 1. Handle Empty Table ---
        if (columns.empty()) {
            std::println("{} \n(Table has no columns)\n", tableName);
            return;
        }

        // --- 2. Determine dimensions and calculate column widths (using visual_width) ---
        size_t num_cols = columns.size();
        size_t num_rows = 0;
        for (const auto& col : columns) {
            num_rows = std::max(num_rows, col->getRowCount());
        }

        std::vector<size_t> col_widths(num_cols);
        for (size_t j = 0; j < num_cols; ++j) {
            // Start with header width
            col_widths[j] = visual_width(columns[j]->getHeader());
            // Find max width in all data cells
            for (size_t i = 0; i < num_rows; ++i) {
                col_widths[j] = std::max(col_widths[j], visual_width(columns[j]->getCellData(i)));
            }
        }

        // --- 3. Print the table using std::print / std::println ---

        // --- Table Title ---
        // NOLINTNEXTLINE(*-fold-init-type)
        const size_t total_width = std::accumulate(col_widths.begin(), col_widths.end(), 0UL) + (num_cols * 3) + 1;
        const size_t title_padding_len = (total_width > visual_width(tableName)) ? (total_width - visual_width(tableName)) / 2 : 0;

        // Print padding, then title
        std::print("{: <{}}", "", title_padding_len); // Left-aligned empty string "" of width title_padding_len
        std::println("{}", tableName);


        // --- Helper to draw horizontal border ---
        auto draw_border = [&]() {
            std::print("+");
            for (const size_t width : col_widths) {
                // std::string(width + 2, '-') is still the easiest way to repeat a char
                std::print("{:-<{}}+", "", width + 2); // Prints '-' repeated (width + 2) times
            }
            std::println("");
        };

        // --- Draw Top Border ---
        draw_border();

        // --- Helper to print a cell with correct padding ---
        auto print_cell = [&](const std::string& text, size_t width) {
            size_t text_width = visual_width(text);
            size_t padding = (width >= text_width) ? (width - text_width) : 0;
            // Print text and then the manual padding
            std::print(" {}{: <{}} |", text, "", padding);
        };

        // --- Draw Header Row ---
        std::print("|");
        for (size_t j = 0; j < num_cols; ++j) {
            print_cell(columns[j]->getHeader(), col_widths[j]);
        }
        std::println("");

        // --- Draw Separator ---
        draw_border();

        // --- Draw Data Rows ---
        for (size_t i = 0; i < num_rows; ++i) {
            std::print("|");
            for (size_t j = 0; j < num_cols; ++j) {
                print_cell(columns[j]->getCellData(i), col_widths[j]);
            }
            std::println("");
        }

        // --- Draw Bottom Border ---
        draw_border();
    }

    inline void print_table(const TableBase& table) {
        print_table(table.getName(), table.getColumns());
    }


    inline void to_csv(const std::string& filename, const std::vector<std::unique_ptr<ColumnBase>>& columns) {
        std::ofstream output(filename);
        if (!output.is_open()) {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        // Write header
        for (size_t j = 0; j < columns.size(); ++j) {
            output << columns[j]->getHeader();
            if (j < columns.size() - 1) {
                output << ",";
            }
        }
        output << "\n";
        // Write data rows
        size_t num_rows = 0;
        for (const auto& col : columns) {
            num_rows = std::max(num_rows, col->getRowCount());
        }

        for (size_t i = 0; i < num_rows; ++i) {
            for (size_t j = 0; j < columns.size(); ++j) {
                output << columns[j]->getCellData(i);
                if (j < columns.size() - 1) {
                    output << ",";
                }
            }
            output << "\n";
        }
        output.close();
    }


    inline nlohmann::json to_json(const std::vector<std::unique_ptr<ColumnBase>>& columns) {
        using json = nlohmann::json;
        json j;
        for (const auto& col : columns) {
            std::vector<std::string> col_data;
            const size_t row_count = col->getRowCount();
            for (size_t i = 0; i < row_count; ++i) {
                col_data.push_back(col->getCellData(i));
            }
            j[col->getHeader()] = col_data;
        }
        return j;
    }

    inline void to_json_file(const std::string& filename, const std::vector<std::vector<std::unique_ptr<ColumnBase>>> &tables, const std::vector<std::string>& tableNames) {
        using json = nlohmann::json;
        json j;
        for (size_t t = 0; t < tables.size(); ++t) {
            j[tableNames[t]] = to_json(tables[t]);
        }
        std::ofstream output(filename);
        if (!output.is_open()) {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        output << j.dump(4);
        output.close();
    }


}