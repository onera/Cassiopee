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
    pointer_message_buffer message_buffer;
    pending_message(pointer_message_buffer pt_msg) : message_buffer(pt_msg) {}
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
  ~PendingMsgContainer() {}

  template <typename... Args>
  void emplace_back(Args const&... args);
  void push_back(pointer_message_buffer& pt_buff);
  void pop(const iterator& it);

  bool empty() const { return m_pending_message.empty(); }
  std::size_t size() const { return m_pending_message.size(); }
  void clear() { m_pending_message.clear(); }

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
