// Sparse matrices in compressed sparse rows (CSR).
//
// The entries of row r are positions [row[r], row[r + 1]) of `col` and of any arrays of
// values the caller keeps beside it (a weight, a delay, the state of a synapse), so that
// one structure carries entries of any kind. compress() builds it from coordinates by a
// stable counting sort, O(rows + entries): entries of a row keep the order in which they
// were given, which makes the result independent of how the coordinates were produced, as
// long as their order is.
//
// Saad 2003, Iterative Methods for Sparse Linear Systems, 2nd ed., §3.4.

#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::linalg {

template <class Index = std::uint32_t>
struct Csr {
  std::size_t rows = 0, cols = 0;
  std::vector<std::uint64_t> row;  // rows + 1 offsets
  std::vector<Index> col;

  std::size_t entries() const { return col.size(); }
  std::uint64_t begin(std::size_t r) const { return row[r]; }
  std::uint64_t end(std::size_t r) const { return row[r + 1]; }
};

// Compresses the entries (r[k], c[k]). If `order` is given it receives the permutation:
// the k-th stored entry was the order[k]-th given, so that values given beside the
// coordinates can be placed by values_out[k] = values_in[order[k]].
template <class Index>
Csr<Index> compress(std::size_t rows, std::size_t cols, std::span<const Index> r, std::span<const Index> c,
                    std::vector<std::uint64_t>* order = nullptr) {
  if (r.size() != c.size()) throw std::invalid_argument("calx::linalg::compress: rows and columns differ in length");
  Csr<Index> m;
  m.rows = rows;
  m.cols = cols;
  m.row.assign(rows + 1, 0);
  for (const Index i : r) {
    if (static_cast<std::size_t>(i) >= rows) throw std::out_of_range("calx::linalg::compress: row out of range");
    m.row[static_cast<std::size_t>(i) + 1]++;
  }
  for (std::size_t i = 0; i < rows; i++) m.row[i + 1] += m.row[i];
  std::vector<std::uint64_t> next(m.row.begin(), m.row.end() - 1);
  m.col.resize(r.size());
  if (order) order->resize(r.size());
  for (std::size_t k = 0; k < r.size(); k++) {
    const std::uint64_t at = next[static_cast<std::size_t>(r[k])]++;
    if (static_cast<std::size_t>(c[k]) >= cols) throw std::out_of_range("calx::linalg::compress: column out of range");
    m.col[at] = c[k];
    if (order) (*order)[at] = k;
  }
  return m;
}

// y = A x for values `a` stored beside the columns.
template <class Index, class Value>
void multiply(const Csr<Index>& A, std::span<const Value> a, std::span<const double> x, std::span<double> y) {
  for (std::size_t i = 0; i < A.rows; i++) {
    double s = 0;
    for (std::uint64_t k = A.row[i]; k < A.row[i + 1]; k++) s += static_cast<double>(a[k]) * x[A.col[k]];
    y[i] = s;
  }
}

}  // namespace calx::linalg
