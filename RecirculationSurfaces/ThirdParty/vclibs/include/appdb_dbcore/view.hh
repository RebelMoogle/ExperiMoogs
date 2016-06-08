//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_VIEW_HH
#define VC_APPDB_VIEW_HH

#include <vector>
#include <unordered_map>

#include "appdb.hh"

//== INCLUDES =================================================================



//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

/** A sequence of nodes.
    \ingroup vc_appdb
 */
class View {
public:
  typedef std::vector<Node*> node_list_t;

  View();
  View(const View& _other);
  virtual ~View();

  virtual void clear();
  virtual void push(Node* _node);

  AbstractDatabase* db() const {
    return !m_nodes.empty() ? m_nodes.front()->db() : 0;
  }
  const node_list_t& nodes() const { return m_nodes; }
  bool includes(Node* _node) const {
    return m_map.find(_node)!=m_map.end();
  }
  int index(Node* _node) const {
    auto ii=m_map.find( _node);
    return ii!=m_map.end() ? int(ii->second) : -1;
  }

  const boost_signals::connection& node_changes() const {
    return m_nc; // e.g., shared_connection_block
  }


  /** @name Utilities
      @{
  */

  /*

    TODO: revert change!
    reformulate to void default template parameter (MSC fails)

  struct  True {
    bool operator()(const Node* _node) const { return true; }
  };

  */

  template <typename ConstIterator,typename Predicate/*=True*/>
  void push_range(const ConstIterator& _begin,const ConstIterator& _end,
                  Predicate _pred/*=True()*/) {
    for (ConstIterator ii=_begin;ii!=_end;++ii)
      if (_pred(*ii))
        push(*ii);
  }
  template <typename ConstIterator>
  void push_range(const ConstIterator& _begin,const ConstIterator& _end) {
    for (ConstIterator ii=_begin;ii!=_end;++ii)
      push(*ii);
  }

  template <typename Predicate/*=True*/>
  void push_children(Node* _node,Predicate& _pred/*=True()*/) {
    std::vector<Node*> ch=_node->children();
    push_range(ch.begin(),ch.end(),_pred);
  }  
  void push_children(Node* _node) {
    std::vector<Node*> ch=_node->children();
    push_range(ch.begin(),ch.end());
  }

  /// @}


protected:
  typedef std::unordered_map<Node*,size_t> node_map_t;

  void connect();
  void disconnect();
  void handle_node_changed(Database::Event _event,Node* _node,void* _data);

  virtual void node_changed(Database::Event _event,Node* _node,void* _data);

  node_list_t               m_nodes;
  node_map_t                m_map;
  boost_signals::connection m_nc;
  int                       m_changing_node;
};

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------


//=============================================================================
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_VIEW_HH defined
