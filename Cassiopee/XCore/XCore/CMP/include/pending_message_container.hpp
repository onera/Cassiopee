/*    
    Copyright 2013-2026 ONERA.

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
// pending_message_container.hpp
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef _CMP_PENDING_MESSAGE_CONTAINER_HPP_
#define _CMP_PENDING_MESSAGE_CONTAINER_HPP_
#include <memory>
#include <vector>
using std::shared_ptr;

namespace CMP {
template <typename MsgBuffer>
class PendingMsgContainer {
 public:
  typedef shared_ptr<MsgBuffer> pointer_message_buffer;
  struct pending_message {
    pointer_message_buffer message_buffer = nullptr;
    pending_message() = default;
    pending_message(pointer_message_buffer pt_msg) : message_buffer(pt_msg) {}
    pending_message(const pending_message& ) = default;
    MsgBuffer& get_message_buffer() { return *message_buffer; }
    const MsgBuffer& get_message_buffer() const { return *message_buffer; }
  };
  typedef std::vector<pending_message> container;
  typedef typename container::iterator iterator;
  typedef typename container::const_iterator const_iterator;

  PendingMsgContainer() : m_pending_message() {}
  PendingMsgContainer(std::size_t maxSize) : m_pending_message() {
    m_pending_message.reserve(maxSize);
  }
  PendingMsgContainer(const PendingMsgContainer&) = default;
  PendingMsgContainer(      PendingMsgContainer&&) = default;
  ~PendingMsgContainer() {}

  PendingMsgContainer& operator = ( const PendingMsgContainer& ) = default;
  PendingMsgContainer& operator = ( PendingMsgContainer     && ) = default;

  template <typename... Args>
  void emplace_back(Args const&... args);
  void push_back(pointer_message_buffer& pt_buff);
  void pop(const iterator& it);

  void resize( std::size_t n ) { m_pending_message.resize(n); }
  bool empty() const { return m_pending_message.empty(); }
  std::size_t size() const { return m_pending_message.size(); }
  void clear() { for ( auto& pm : m_pending_message ) pm.message_buffer->clear(); }

  void waitAll();
  // If iterator is equal to end(), not complete message found
  iterator get_first_complete_message();

  pending_message& operator[](std::size_t i);
  const pending_message& operator[](std::size_t i) const;

  pending_message& back() { return m_pending_message.back(); }
  const pending_message& back() const { return m_pending_message.back(); }
  MsgBuffer& back_message_buffer() {
    return *m_pending_message.back().message_buffer;
  }
  const MsgBuffer& back_message_buffer() const {
    return *m_pending_message.back().message_buffer;
  }

  pending_message& front() { return m_pending_message.front(); }
  const pending_message& front() const { return m_pending_message.front(); }
  MsgBuffer& front_message_buffer() {
    return *m_pending_message.front().message_buffer;
  }
  const MsgBuffer& front_message_buffer() const {
    return *m_pending_message.front().message_buffer;
  }

  iterator begin() { return m_pending_message.begin(); }
  const_iterator begin() const { return cbegin(); }
  const_iterator cbegin() const { return m_pending_message.begin(); }
  iterator end() { return m_pending_message.end(); }
  const_iterator end() const { return cend(); }
  const_iterator cend() const { return m_pending_message.end(); }

 private:
  container m_pending_message;
};
}  // namespace CMP
#endif
