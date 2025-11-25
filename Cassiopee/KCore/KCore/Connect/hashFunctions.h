/*    
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstddef>
#include <utility>

/* Hash functions */
// Template parameter is a topology struct (see topologyMapping.h), unless
// otherwise stated
// BernsteinHash: Bernstein djb2 hash function
template <typename T>
struct BernsteinHash
{
  std::size_t operator()(const T& key) const
  {
    std::size_t hash = 5381;
    for (std::size_t i = 0; i < key.n_; i++)
    {
      hash = ((hash << 5) + hash) + (std::size_t)key.p_[i];
    }
    return hash;
  }
};

// FowlerNollVoHash: FNV-1a (Fowler-Noll-Vo) hash function where the prime
// multiplier is the golden ratio (32-bit)
template <typename T>
struct FowlerNollVoHash
{
  std::size_t operator()(const T& key) const
  {
    std::size_t hash = 0x811c9dc5; // 32-bit offset basis
    for (std::size_t i = 0; i < key.n_; i++)
    {
      hash = ((std::size_t)key.p_[i] ^ hash) * 0x9e3779b9;
    }
    return hash;
  }
};

// Hash function for std::pair of integers based on FNV-1a
struct pairHash
{
  std::size_t operator () (const std::pair<E_Int, E_Int>& pair) const
  {
    std::size_t hash = 0x811c9dc5; // 32-bit offset basis
    hash = (pair.first ^ hash) * 0x9e3779b9;
    hash = (pair.second ^ hash) * 0x9e3779b9;
    return hash;
  }
};

// JenkinsHash: Jenkins's one-at-a-time hash function
template <typename T>
struct JenkinsHash
{
  std::size_t operator()(const T& key) const
  {
    std::size_t hash = 0;
    for (std::size_t i = 0; i < key.n_; i++)
    {
      hash += (std::size_t)key.p_[i];
      hash += hash << 10;
      hash ^= hash >> 6;
    }
    hash += hash << 3;
    hash ^= hash >> 11;
    hash += hash << 15;
    return hash;
  }
};
