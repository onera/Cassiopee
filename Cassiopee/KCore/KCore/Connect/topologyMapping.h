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
#include "kcore.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>

/* Topology structs */
// Topology object designed to work with std::map and std::unordered_map
// No explicit limit on the number of items contained in the Topology
struct Topology
{
  E_Bool isDegen_; // whether the Topology is degenerated
  std::size_t n_; // number of items contained in the Topology
  std::size_t size_; // number of unique items contained in the Topology
  std::vector<E_Int> p_;

  Topology() {}
  Topology(const E_Int* p, std::size_t n, E_Bool search4Degen=false)
  {
    set(p, n, search4Degen);
  }
  Topology(const std::vector<E_Int>& p, E_Bool search4Degen=false)
  {
    set(p, search4Degen);
  }
  
  void set(const E_Int* p, std::size_t n, E_Bool search4Degen=false)
  {
    size_ = 0; isDegen_ = false;
    n_ = n; p_.assign(p, p+n_);
    std::sort(p_.begin(), p_.end());
    if (search4Degen) countUnique();
  }
  
  void set(const std::vector<E_Int>& p, E_Bool search4Degen=false)
  {
    size_ = 0; isDegen_ = false;
    n_ = p.size(); p_ = p;
    std::sort(p_.begin(), p_.end());
    if (search4Degen) countUnique();
  }

  void countUnique()
  {
    size_ = 0; isDegen_ = false;
    std::unordered_set<E_Int> pSet;
    for (std::size_t i = 0; i < n_; i++)
    {
      if (pSet.insert(p_[i]).second) size_++;
    }
    if (size_ != n_) isDegen_ = true;
  }

  E_Bool operator<(const Topology& other) const
  {
    std::size_t minsize = std::min(n_, other.n_);
    for (std::size_t i = 0; i < minsize; i++)
    {
      if (p_[i] < other.p_[i]) return true;
      else if (p_[i] > other.p_[i]) return false;
    }

    // If all common items are equal, the shorter vector is considered
    // smaller
    return n_ < other.n_;
  }

  E_Bool operator==(const Topology& other) const
  {
    if (n_ != other.n_) return false;
    for (std::size_t i = 0; i < n_; i++)
    {
      if (p_[i] != other.p_[i]) return false;
    }
    return true;
  }

  friend std::ostream& operator<<(std::ostream& os, const Topology& T)
  {
    os << "Topology: \n";
    os << "  - number of items: " << T.n_ << "\n";
    os << "  - number of unique items: " << T.size_ << "\n";
    os << "  - is degenerated?: " << T.isDegen_ << "\n";
    os << "  - data: ";
    for (std::size_t i = 0; i < T.n_; i++) os << T.p_[i] << " ";
    os << "\n";
    return os;
  }
};

// Same as Topology but optimised (no dynamic allocation)
const std::size_t nmaxitems = 9; // nbre d'index par item
struct TopologyOpt
{
  E_Bool isDegen_; // whether the Topology is degenerated
  std::size_t n_; // number of items contained in the Topology
  std::size_t size_; // number of unique items contained in the Topology
  E_Int p_[nmaxitems]; // indices

  TopologyOpt() { for (size_t i = 0; i < nmaxitems; i++) p_[i] = -1; }
  TopologyOpt(const E_Int* p, const std::size_t n, E_Bool search4Degen=false)
  {
    set(p, n, search4Degen);
  }
  TopologyOpt(const std::vector<E_Int>& p, const std::size_t n,
              E_Bool search4Degen=false)
  {
    set(p, n, search4Degen);
  }
  
  void set(const E_Int* p, const std::size_t n, E_Bool search4Degen=false)
  {
    size_ = 0; isDegen_ = false;
    n_ = n; assert(n_ <= nmaxitems);
    //for (std::size_t i = 0; i < n; i++) p_[i] = p[i];
    memcpy(p_, p, n_*sizeof(E_Int));
    std::sort(p_, p_ + n_);
    if (search4Degen) countUnique();
  }
  
  void set(const std::vector<E_Int>& p, const std::size_t n,
           E_Bool search4Degen=false)
  {
    size_ = 0; isDegen_ = false;
    n_ = n; assert(n_ <= nmaxitems);
    for (std::size_t i = 0; i < n_; i++) p_[i] = p[i];
    std::sort(p_, p_ + n_);
    if (search4Degen) countUnique();
  }

  void countUnique()
  {
    size_ = 0; isDegen_ = false;
    std::unordered_set<E_Int> pSet;
    for (std::size_t i = 0; i < n_; i++)
    {
      if (pSet.insert(p_[i]).second) size_++;
    }
    if (size_ != n_) isDegen_ = true;
  }

  E_Bool operator<(const TopologyOpt& other) const
  {
    std::size_t minsize = std::min(n_, other.n_);
    for (std::size_t i = 0; i < minsize; i++)
    {
      if (p_[i] < other.p_[i]) return true;
      else if (p_[i] > other.p_[i]) return false;
    }

    // If all common items are equal, the shorter array is considered smaller
    return n_ < other.n_;
  }

  E_Bool operator==(const TopologyOpt& other) const
  {
    if (n_ != other.n_) return false;
    for (std::size_t i = 0; i < n_; i++)
    {
      if (p_[i] != other.p_[i]) return false;
    }
    return true;
  }

  friend std::ostream& operator<<(std::ostream& os, const TopologyOpt& T)
  {
    os << "TopologyOpt: \n";
    os << "  - number of items: " << T.n_ << "\n";
    os << "  - number of unique items: " << T.size_ << "\n";
    os << "  - is degenerated?: " << T.isDegen_ << "\n";
    os << "  - data: ";
    for (std::size_t i = 0; i < T.n_; i++) os << T.p_[i] << " ";
    os << "\n";
    return os;
  }
};
