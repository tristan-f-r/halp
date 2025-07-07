"""
.. module:: mixed_hypergraph
   :synopsis: Defines MixedHypergraph class for the basic properties
            of a mixed hypergraph, along with the relevant structures
            regarding nodes, hyperedges, adjacency, etc.

"""

import copy

from halp.undirected_hypergraph import UndirectedHypergraph
from halp.directed_hypergraph import DirectedHypergraph

class MixedHypergraph(object):
    """
    The MixedHypergraph class provides a mixed hypergraph object
    and associated functions for basic properties of mixed hypergraphs.

    A mixed hypergraph contains nodes and hyperedges. Each hyperedge
    can either be a 'directed hyperedge' which connects a tail set of
    nodes to a head set of nodes (where the tail and head cannot both be empty),
    or an undirected hyperedge, which connects a set of nodes.

    A node is simply any hashable type. See "add_node" or "add_nodes" for
    more details.

    A directed hyperedge is a tuple of the tail nodes and the head nodes.
    An undirected hyperedge is any iterable container of nodes.
    This class assigns (upon adding) and refers to each hyperedge by an
    internal ID. See "add_hyperedge" or "add_hyperedges" for more details.

    Directed self-loops are allowed, but parallel (multi) hyperedges are not.

    :note: This class uses several data structures to store a mixed
        hypergraph. Since these structures must stay in sync (see: __init__),
        we highly recommend that only the public methods be used for accessing
        and modifying the hypergraph.

    Examples:
    Create an empty mixed hypergraph (no nodes or hyperedges):

    >>> H = MixedHypergraph()

    Add nodes (with or without attributes) to the hypergraph
    one at a time (see: add_node) or several at a time (see: add_nodes):

    >>> H.add_nodes(["A", "B", "C", "D"], {color: "black"})
    
    Add hyperedges (with or without attributes) to the hypergraph one
    at a time (see: add_hyperedge) or several at a time (see: add_hyperedges):

    >>> H.add_hyperedges((["A"], ["B"]), (["A", "B"], ["C", "D"]))

    Update attributes of existing nodes and hyperedges by simulating adding the
    node or hyperedge again, with the [potentially new] attribute appended:

    >>> H.add_node("A", label="sink")
    >>> H.add_hyperedge((["A", "B"], ["C", "D"]), weight=5)

    """

    def __init__(self):
        """Constructor for the MixedHypergraph class.
        Initializes all internal data structures used for the rapid
        execution of most of the fundamental hypergraph queries.

        """
        # _node_attributes: a dictionary mapping a node (any hashable type)
        # to a dictionary of attributes of that node.
        #
        # Provides O(1) time access to the attributes of a node.
        #
        # Used in the implementation of methods such as add_node and
        # get_node_attributes.
        #
        self._node_attributes = {}

        # _hyperedge_attributes: a dictionary mapping a hyperedge ID
        # (initially created by the call to add_hyperedge or add_hyperedges)
        # to a dictionary of attributes of that hyperedge.
        # Given a hyperedge ID, _hyperedge_attributes[hyperedge_id] stores
        # the tail of the hyperedge as specified by the user (as "tail"),
        # the head of the hyperedge as specified by the user (as "head"),
        # and the weight of the hyperedge (as "weight").
        # For internal purposes, it also stores the frozenset versions of
        # the tail and head (as "__frozen_tail" and "__frozen_head").
        #
        # Provides O(1) time access to the attributes of a hyperedge.
        #
        # Used in the implementation of methods such as add_hyperedge and
        # get_hyperedge_attributes.
        #
        self._hyperedge_attributes = {}

        # The star of a node is the set of undirected hyperedges the
        # node is a member of.
        #
        # _star: a dictionary mapping a node to the set of undirected
        # hyperedges that are in that node's star.
        #
        # Provides O(1) time access to a reference to the set of hyperedges
        # that a node is a member of.
        #
        # Used in the implementation of methods such as add_node and
        # remove_hyperedge.
        #
        self._star = {}

        # The forward star of a node is the set of directed hyperedges
        # such that the node is in the tail of each hyperedge in that set.
        #
        # _forward_star: a dictionary mapping a node to the set of hyperedges
        # that are in that node's forward star.
        #
        # Provides O(1) time access to a reference to the set of outgoing
        # hyperedges from a node.
        #
        # Used in the implementation of methods such as add_node and
        # remove_hyperedge.
        #
        self._forward_star = {}

        # The backward star of a node is the set of directed hyperedges such
        # that the node is in the head of each hyperedge in that set.
        #
        # _backward_star: a dictionary mapping a node to the set of hyperedges
        # that are in that node's backward star.
        #
        # Provides O(1) time access to a reference to the set of incoming
        # hyperedges from a node.
        #
        # Used in the implementation of methods such as add_node and
        # remove_hyperedge.
        #
        self._backward_star = {}

        # _successors: a 2-dimensional dictionary mapping (first) a tail set
        # and (second) a head set of a directed hyperedge to the ID of the
        # corresponding hyperedge. We represent each tail set and each head set
        # by a frozenset, so that the structure is hashable.
        #
        # Provides O(1) time access to the ID of the of the hyperedge
        # connecting a specific tail frozenset to a specific head frozenset.
        # Given a tail frozenset, it also provides O(1) time access to a
        # reference to the dictionary mapping head frozensets to hyperedge IDs;
        # these hyperedges are precisely those in the forward star of this
        # tail frozenset.
        #
        self._successors = {}

        # _predecessors: a 2-dimensional dictionary mapping (first) a head set
        # and (second) a tail set of a directed hyperedge to the ID of the
        # corresponding hyperedge. We represent each tail set and each head set by a
        # frozenset, so that the structure is hashable.
        #
        # Provides O(1) time access to the ID of the of the hyperedge
        # connecting a specific head frozenset to a specific tail frozenset.
        # Given a head frozenset, it also provides O(1) time access to a
        # reference to the dictionary mapping tail frozensets to hyperedge IDs;
        # these hyperedges are precisely those in the backward star of this
        # head frozenset.
        #
        self._predecessors = {}

        # _node_set_to_hyperedge: a dictionary mapping a set of nodes to the ID
        # of the undirected hyperedge they compose. We represent the node set by a
        # frozenset, so that the structure is hashable.
        #
        # Provides O(1) time access to the ID of the hyperedge that
        # a specific frozenset of nodes composes.
        #
        self._node_set_to_hyperedge = {}

        # _current_hyperedge_id: an int representing the hyperedge ID that
        # was most recently assigned by the class (since users don't
        # name/ID their own hyperedges); hyperedges being added are issued
        # ID "e"+_current_hyperedge_id.
        #
        # Since the class takes responsibility for giving hyperedges
        # their IDs (i.e. a unique identifier; could be alternatively viewed
        # as a unique name, label, etc.), the issued IDs need to be kept
        # track of. A consecutive issuing of integer IDs to the hyperedges is a
        # simple strategy to ensure their uniqueness and allow for
        # intuitive readability.
        #
        # e.g., _current_hyperedge_id = 4  implies that 4 hyperedges have
        # been added to the hypergraph, and that "e4" was the most recently
        # assigned hyperedge.
        #
        # Note: An hyperedge, once added, will receive a unique ID. If this
        # hyperedge is removed and subsequently re-added, it will not receive
        # the same ID as it was issued when it was originally added.
        #
        self._current_hyperedge_id = 0

    def _combine_attribute_arguments(self, attr_dict, attr):
        """Combines attr_dict and attr dictionaries, by updating attr_dict
            with attr.

        :param attr_dict: dictionary of attributes of the node.
        :param attr: keyword arguments of attributes of the node;
                    attr's values will override attr_dict's values
                    if both are provided.
        :returns: dict -- single dictionary of [combined] attributes.
        :raises: AttributeError -- attr_dict argument must be a dictionary.

        """
        # If no attribute dict was passed, treat the keyword
        # arguments as the dict
        if attr_dict is None:
            attr_dict = attr
        # Otherwise, combine the passed attribute dict with
        # the keyword arguments
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise AttributeError("attr_dict argument \
                                     must be a dictionary.")
        return attr_dict

    def has_node(self, node):
        """Determines if a specific node is present in the hypergraph.

        :param node: reference to the node whose presence is being checked.
        :returns: bool -- true iff the node exists in the hypergraph.

        """
        return node in self._node_attributes

    def add_node(self, node, attr_dict=None, **attr):
        """Adds a node to the graph, along with any related attributes
           of the node.

        :param node: reference to the node being added.
        :param attr_dict: dictionary of attributes of the node.
        :param attr: keyword arguments of attributes of the node;
                    attr's values will override attr_dict's values
                    if both are provided.
        :note: Following the NetworkX model, we allow both a dictionary
            of attributes to be passed in by the user as well as key-value
            pairs of attributes to be passed in by the user to provide
            greater flexibility. This pattern is followed in methods such
            as add_nodes, add_hyperedge, and add_hyperedges.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> attributes = {label: "positive"}
            >>> H.add_node("A", attributes)
            >>> H.add_node("B", label="negative")
            >>> H.add_node("C", attributes, root=True)

        """
        attr_dict = self._combine_attribute_arguments(attr_dict, attr)

        # If the node hasn't previously been added, add it along
        # with its attributes
        if not self.has_node(node):
            self._node_attributes[node] = attr_dict
            self._forward_star[node] = set()
            self._backward_star[node] = set()
            self._star[node] = set()
        # Otherwise, just update the node's attributes
        else:
            self._node_attributes[node].update(attr_dict)

    def add_nodes(self, nodes, attr_dict=None, **attr):
        """Adds multiple nodes to the graph, along with any related attributes
            of the nodes.

        :param nodes: iterable container to either references of the nodes
                    OR tuples of (node reference, attribute dictionary);
                    if an attribute dictionary is provided in the tuple,
                    its values will override both attr_dict's and attr's
                    values.
        :param attr_dict: dictionary of attributes shared by all the nodes.
        :param attr: keyword arguments of attributes of the node;
                    attr's values will override attr_dict's values
                    if both are provided.

        See also:
        add_node

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> attributes = {label: "positive"}
            >>> node_list = ["A",
                             ("B", {label="negative"}),
                             ("C", {root=True})]
            >>> H.add_nodes(node_list, attributes)

        """
        attr_dict = self._combine_attribute_arguments(attr_dict, attr)

        for node in nodes:
            # Note: This won't behave properly if the node is actually a tuple
            if type(node) is tuple:
                # See ("B", {label="negative"}) in the documentation example
                new_node, node_attr_dict = node
                # Create a new dictionary and load it with node_attr_dict and
                # attr_dict, with the former (node_attr_dict) taking precedence
                new_dict = attr_dict.copy()
                new_dict.update(node_attr_dict)
                self.add_node(new_node, new_dict)
            else:
                # See "A" in the documentation example
                self.add_node(node, attr_dict.copy())

    def remove_node(self, node):
        """Removes a node and its attributes from the hypergraph. Removes
        every hyperedge that contains this node in either the head or the tail.

        :param node: reference to the node being added.
        :raises: ValueError -- No such node exists.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> H.add_node("A", label="positive")
            >>> H.remove_node("A")

        """
        if not self.has_node(node):
            raise ValueError("No such node exists.")

        # Loop over every undirected hyperedge in the star of the node;
        # i.e., over every undirected hyperedge that contains the node
        for hyperedge_id in self._star[node]:
            frozen_nodes = \
                self._hyperedge_attributes[hyperedge_id]["__frozen_nodes"]
            # Remove the node set composing the hyperedge
            del self._node_set_to_hyperedge[frozen_nodes]
            # Remove this hyperedge's attributes
            del self._hyperedge_attributes[hyperedge_id]

        # Remove every directed hyperedge which is in the forward star of the node
        forward_star = self.get_forward_star(node)
        for hyperedge_id in forward_star:
            self.remove_hyperedge(hyperedge_id)

        # Remove every directed  hyperedge which is in the backward star of the node
        # but that is not also in the forward start of the node (to handle
        # overlapping hyperedges)
        backward_star = self.get_backward_star(node)
        for hyperedge_id in backward_star - forward_star:
            self.remove_hyperedge(hyperedge_id)

        # Remove all of the node's stars
        del self._forward_star[node]
        del self._backward_star[node]
        del self._star[node]

        # Remove node's attributes dictionary
        del self._node_attributes[node]

    def remove_nodes(self, nodes):
        """Removes multiple nodes and their attributes from the graph. If
        the nodes are part of any hyperedges, those hyperedges are removed
        as well.

        :param nodes: iterable container to either references of the nodes
                    OR tuples of (node reference, attribute dictionary);
                    if an attribute dictionary is provided in the tuple,
                    its values will override both attr_dict's and attr's
                    values.

        See also:
        remove_node

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> attributes = {label: "positive"}
            >>> node_list = ["A",
                             ("B", {label="negative"}),
                             ("C", {root=True})]
            >>> H.add_nodes(node_list, attributes)
            >>> H.remove_nodes(["A", "B", "C"])

        """
        for node in nodes:
            self.remove_node(node)

    def trim_node(self, node):
        """Removes a node from the hypergraph. Modifies hyperedges with the 
        trimmed node in their head or tail so that they no longer include 
        the trimmed node. If a hyperedge has solely the trimmed node in its
        head or tail, that hyperedge is removed.
        
        Note: hyperedges modified this way will have different IDs than before
        
        Examples:
        ::
        
            >>> H = MixedHypergraph()
            >>> node_list = ['A', 'B', 'C', 'D']
            >>> H.add_nodes(node_list)
            >>> H.add_hyperedge(['A','B'],['C','D'])
            >>> H.trim_node('A')
        """
    
        def get_attrs(H, hyperedge):
            #copies the attribute dictionary of a hyperedge except for the head, tail, and nodes
            new_attrs = {}
            old_attrs = H.get_hyperedge_attributes(hyperedge)
        
            for key in old_attrs:
                if key not in {'head', 'tail', 'nodes'}:
                    new_attrs[key] = old_attrs[key]
            return new_attrs

        remove_set = set()
        # Handle undirected hyperedges
        s = self.get_star(node)
    
        for hedge in s:
            nodes = set(self.get_undirected_hyperedge_nodes(hedge))
            if len(nodes) > 1:
                new_nodes = nodes - {node}
                attrs = get_attrs(self, hedge)
                self.add_undirected_hyperedge(new_nodes, attrs)
            remove_set.add(hedge)
        
        # Handle directed hyperedges
        fs = self.get_forward_star(node)
        bs = self.get_backward_star(node)
    
        for hedge in fs:
            tail = set(self.get_directed_hyperedge_tail(hedge))
            head = set(self.get_directed_hyperedge_head(hedge))
            if len(tail) > 1:
                new_tail = tail - {node}
                attrs = get_attrs(self, hedge)
                self.add_directed_hyperedge(new_tail, head, attrs)
            remove_set.add(hedge)
            
        for hedge in bs:
            tail = set(self.get_directed_hyperedge_tail(hedge))
            head = set(self.get_directed_hyperedge_head(hedge))
            if len(head) > 1:
                new_head = head - {node}
                attrs = get_attrs(self, hedge)
                self.add_directed_hyperedge(tail, new_head, attrs)
            remove_set.add(hedge)

        for hedge in remove_set:
            self.remove_hyperedge(hedge)
            
        self.remove_node(node)
    
    def trim_nodes(self, nodes):
        """Trims multiple nodes from the hypergraph (see trim_node() for details)"""
        for node in nodes:
            self.trim_node(node)
            
    def get_node_set(self):
        """Returns the set of nodes that are currently in the hypergraph.

        :returns: set -- all nodes currently in the hypergraph

        """
        return set(self._node_attributes.keys())

    def node_iterator(self):
        """Provides an iterator over the nodes.

        """
        return iter(self._node_attributes)

    def get_node_attribute(self, node, attribute_name):
        """Given a node and the name of an attribute, get a copy
        of that node's attribute.

        :param node: reference to the node to retrieve the attribute of.
        :param attribute_name: name of the attribute to retrieve.
        :returns: attribute value of the attribute_name key for the
                specified node.
        :raises: ValueError -- No such node exists.
        :raises: ValueError -- No such attribute exists.

        """
        if not self.has_node(node):
            raise ValueError("No such node exists.")
        elif attribute_name not in self._node_attributes[node]:
            raise ValueError("No such attribute exists.")
        else:
            return copy.\
                copy(self._node_attributes[node][attribute_name])

    def get_node_attributes(self, node):
        """Given a node, get a dictionary with copies of that node's
        attributes.

        :param node: reference to the node to retrieve the attributes of.
        :returns: dict -- copy of each attribute of the specified node.
        :raises: ValueError -- No such node exists.

        """
        if not self.has_node(node):
            raise ValueError("No such node exists.")
        attributes = {}
        for attr_name, attr_value in self._node_attributes[node].items():
            attributes[attr_name] = copy.copy(attr_value)
        return attributes

    def _assign_next_hyperedge_id(self):
        """Returns the next [consecutive] ID to be assigned
            to a hyperedge.
        :returns: str -- hyperedge ID to be assigned.

        """
        self._current_hyperedge_id += 1
        return "e" + str(self._current_hyperedge_id)

    def add_undirected_hyperedge(self, nodes, attr_dict=None, **attr):
        """Adds an undirected hyperedge to the hypergraph, along with any
            related attributes of the hyperedge.
            This method will automatically add any node from the node set
            that was not in the hypergraph.
            A hyperedge without a "weight" attribute specified will be
            assigned the default value of 1.

        :param nodes: iterable container of references to nodes in the
                    hyperedge to be added.
        :param attr_dict: dictionary of attributes of the hyperedge being
                        added.
        :param attr: keyword arguments of attributes of the hyperedge;
                    attr's values will override attr_dict's values
                    if both are provided.
        :returns: str -- the ID of the hyperedge that was added.
        :raises: ValueError -- nodes arguments cannot be empty.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> x = H.add_undirected_hyperedge(["A", "B", "C"])
            >>> y = H.add_undirected_hyperedge(("A", "D"), weight=2)
            >>> z = H.add_undirected_hyperedge(set(["B", "D"]), {color: "red"})

        """
        attr_dict = self._combine_attribute_arguments(attr_dict, attr)

        # Don't allow empty node set (invalid hyperedge)
        if not nodes:
            raise ValueError("nodes argument cannot be empty.")

        # Use frozensets for node sets to allow for hashable keys
        frozen_nodes = frozenset(nodes)

        is_new_hyperedge = not self.has_undirected_hyperedge(frozen_nodes)
        if is_new_hyperedge:
            # Add nodes to graph (if not already present)
            self.add_nodes(frozen_nodes)

            # Create new hyperedge name to use as reference for that hyperedge
            hyperedge_id = self._assign_next_hyperedge_id()

            # For each node in the node set, add hyperedge to the node's star
            for node in frozen_nodes:
                self._star[node].add(hyperedge_id)

            # Add the hyperedge ID as the hyperedge that the node set composes
            self._node_set_to_hyperedge[frozen_nodes] = hyperedge_id

            # Assign some special attributes to this hyperedge. We assign
            # a default weight of 1 to the hyperedge. We also store the
            # original node set in order to return them exactly as the
            # user passed them into add_hyperedge.
            self._hyperedge_attributes[hyperedge_id] = \
                {"nodes": nodes, "__frozen_nodes": frozen_nodes, "weight": 1}
        else:
            # If its not a new hyperedge, just get its ID to update attributes
            hyperedge_id = self._node_set_to_hyperedge[frozen_nodes]

        # Set attributes and return hyperedge ID
        self._hyperedge_attributes[hyperedge_id].update(attr_dict)
        return hyperedge_id

    def add_directed_hyperedge(self, tail, head, attr_dict=None, **attr):
        """Adds a directed hyperedge to the hypergraph, along with any
        related attributes of the hyperedge.
        This method will automatically add any node from the tail and
        head that was not in the hypergraph.
        A hyperedge without a "weight" attribute specified will be
        assigned the default value of 1.

        :param tail: iterable container of references to nodes in the
                    tail of the hyperedge to be added.
        :param head: iterable container of references to nodes in the
                    head of the hyperedge to be added.
        :param attr_dict: dictionary of attributes shared by all
                    the hyperedges.
        :param attr: keyword arguments of attributes of the hyperedge;
                    attr's values will override attr_dict's values
                    if both are provided.
        :returns: str -- the ID of the hyperedge that was added.
        :raises: ValueError -- tail and head arguments cannot both be empty.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> x = H.add_directed_hyperedge(["A", "B"], ["C", "D"])
            >>> y = H.add_directed_hyperedge(("A", "C"), ("B"), 'weight'=2)
            >>> z = H.add_directed_hyperedge(set(["D"]),
                                             set(["A", "C"]),
                                             {color: "red"})

        """
        attr_dict = self._combine_attribute_arguments(attr_dict, attr)

        # Don't allow both empty tail and head containers (invalid hyperedge)
        if not tail and not head:
            raise ValueError("tail and head arguments \
                             cannot both be empty.")

        # Use frozensets for tail and head sets to allow for hashable keys
        frozen_tail = frozenset(tail)
        frozen_head = frozenset(head)

        # Initialize a successor dictionary for the tail and head, respectively
        if frozen_tail not in self._successors:
            self._successors[frozen_tail] = {}
        if frozen_head not in self._predecessors:
            self._predecessors[frozen_head] = {}

        is_new_hyperedge = not self.has_directed_hyperedge(frozen_tail, frozen_head)
        if is_new_hyperedge:
            # Add tail and head nodes to graph (if not already present)
            self.add_nodes(frozen_head)
            self.add_nodes(frozen_tail)

            # Create new hyperedge name to use as reference for that hyperedge
            hyperedge_id = self._assign_next_hyperedge_id()

            # Add hyperedge to the forward-star and to the backward-star
            # for each node in the tail and head sets, respectively
            for node in frozen_tail:
                self._forward_star[node].add(hyperedge_id)
            for node in frozen_head:
                self._backward_star[node].add(hyperedge_id)

            # Add the hyperedge as the successors and predecessors
            # of the tail set and head set, respectively
            self._successors[frozen_tail][frozen_head] = hyperedge_id
            self._predecessors[frozen_head][frozen_tail] = hyperedge_id

            # Assign some special attributes to this hyperedge. We assign
            # a default weight of 1 to the hyperedge. We also store the
            # original tail and head sets in order to return them exactly
            # as the user passed them into add_hyperedge.
            self._hyperedge_attributes[hyperedge_id] = \
                {"tail": tail, "__frozen_tail": frozen_tail,
                 "head": head, "__frozen_head": frozen_head,
                 "weight": 1}
        else:
            # If its not a new hyperedge, just get its ID to update attributes
            hyperedge_id = self._successors[frozen_tail][frozen_head]

        # Set attributes and return hyperedge ID
        self._hyperedge_attributes[hyperedge_id].update(attr_dict)
        return hyperedge_id

    def add_undirected_hyperedges(self, hyperedges, attr_dict=None, **attr):
        """Adds multiple undirected hyperedges to the graph, along with any
            related attributes of the hyperedges.
            If any node of a hyperedge has not previously been added to the
            hypergraph, it will automatically be added here.
            Hyperedges without a "weight" attribute specified will be
            assigned the default value of 1.

        :param hyperedges: iterable container to references of the node sets
        :param attr_dict: dictionary of attributes shared by all
                    the hyperedges being added.
        :param attr: keyword arguments of attributes of the hyperedges;
                    attr's values will override attr_dict's values
                    if both are provided.
        :returns: list -- the IDs of the hyperedges added in the order
                    specified by the hyperedges container's iterator.

        See also:
        add_hyperedge

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> hyperedge_list = (["A", "B", "C"],
                                  ("A", "D"),
                                  set(["B", "D"]))
            >>> hyperedge_ids = H.add_hyperedges(hyperedge_list)

        """
        attr_dict = self._combine_attribute_arguments(attr_dict, attr)

        hyperedge_ids = []

        for nodes in hyperedges:
            hyperedge_id = self.add_undirected_hyperedge(nodes, attr_dict.copy())
            hyperedge_ids.append(hyperedge_id)

        return hyperedge_ids

    def add_directed_hyperedges(self, hyperedges, attr_dict=None, **attr):
        """Adds multiple directed hyperedges to the graph, along with any
            related attributes of the hyperedges.
            If any node in the tail or head of any hyperedge has not
            previously been added to the hypergraph, it will automatically
            be added here. Hyperedges without a "weight" attribute specified
            will be assigned the default value of 1.

        :param hyperedges: iterable container to either tuples of
                    (tail reference, head reference) OR tuples of
                    (tail reference, head reference, attribute dictionary);
                    if an attribute dictionary is provided in the tuple,
                    its values will override both attr_dict's and attr's
                    values.
        :param attr_dict: dictionary of attributes shared by all
                    the hyperedges.
        :param attr: keyword arguments of attributes of the hyperedges;
                    attr's values will override attr_dict's values
                    if both are provided.
        :returns: list -- the IDs of the hyperedges added in the order
                    specified by the hyperedges container's iterator.

        See also:
        add_hyperedge

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> xyz = hyperedge_list = ((["A", "B"], ["C", "D"]),
                                        (("A", "C"), ("B"), {'weight': 2}),
                                        (set(["D"]), set(["A", "C"])))
            >>> H.add_directed_hyperedges(hyperedge_list)

        """
        attr_dict = self._combine_attribute_arguments(attr_dict, attr)

        hyperedge_ids = []

        for hyperedge in hyperedges:
            if len(hyperedge) == 3:
                # See ("A", "C"), ("B"), {weight: 2}) in the
                # documentation example
                tail, head, hyperedge_attr_dict = hyperedge
                # Create a new dictionary and load it with node_attr_dict and
                # attr_dict, with the former (node_attr_dict) taking precedence
                new_dict = attr_dict.copy()
                new_dict.update(hyperedge_attr_dict)
                hyperedge_id = self.add_directed_hyperedge(tail, head, new_dict)
            else:
                # See (["A", "B"], ["C", "D"]) in the documentation example
                tail, head = hyperedge
                hyperedge_id = \
                    self.add_directed_hyperedge(tail, head, attr_dict.copy())
            hyperedge_ids.append(hyperedge_id)

        return hyperedge_ids

    def is_hyperedge_id_undirected(self, hyperedge_id):
        """Determines if a hyperedge id is directed.
        
        :param hyperedge_id: ID of the hyperedge whose directionality is
                            being checked.
        
        :returns: bool -- true iff the hyperedge whose id is hyperedge_id is directed."""
        return self.has_hyperedge_id(hyperedge_id) and \
            "__frozen_nodes" in self._hyperedge_attributes[hyperedge_id]


    def is_hyperedge_id_directed(self, hyperedge_id):
        """Determines if a hyperedge id is directed.
        
        :param hyperedge_id: ID of the hyperedge whose directionality is
                            being checked.
        
        :returns: bool -- true iff the hyperedge whose id is hyperedge_id is directed."""
        return self.has_hyperedge_id(hyperedge_id) and \
            not self.is_hyperedge_id_undirected(hyperedge_id)

    def remove_hyperedge(self, hyperedge_id):
        """Removes a hyperedge and its attributes from the hypergraph.

        :param hyperedge_id: ID of the hyperedge to be removed.
        :raises: ValueError -- No such hyperedge exists.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> xyz = hyperedge_list = ((["A"], ["B", "C"]),
                                        (("A", "B"), ("C"), {'weight': 2}),
                                        (set(["B"]), set(["A", "C"])))
            >>> H.add_directed_hyperedges(hyperedge_list)
            >>> H.remove_hyperedge(xyz[0])
            >>> x = H.add_undirected_hyperedge(["A", "B", "C"])
            >>> H.remove_hyperedge(x)

        """
        if not self.has_hyperedge_id(hyperedge_id):
            raise ValueError("No such hyperedge exists.")

        hyperedge_attribute = self._hyperedge_attributes[hyperedge_id]

        if self.is_hyperedge_id_undirected(hyperedge_id):
            frozen_nodes = hyperedge_attribute["__frozen_nodes"]
            # We have an undirected hyperedge
            # Remove this hyperedge from the star of every node in the hyperedge
            for node in frozen_nodes:
                self._star[node].remove(hyperedge_id)

            # Remove this set as the composer of the hyperedge
            del self._node_set_to_hyperedge[frozen_nodes]
        else:
            # We have a directed hyperedge
            frozen_tail = hyperedge_attribute["__frozen_tail"]
            frozen_head = hyperedge_attribute["__frozen_head"]

            # Remove this hyperedge from the forward-star of every tail node
            for node in frozen_tail:
                self._forward_star[node].remove(hyperedge_id)
            # Remove this hyperedge from the backward-star of every head node
            for node in frozen_head:
                self._backward_star[node].remove(hyperedge_id)

            # Remove frozen_head as a successor of frozen_tail
            del self._successors[frozen_tail][frozen_head]
            # If that tail is no longer the tail of any hyperedge, remove it
            # from the successors dictionary
            if self._successors[frozen_tail] == {}:
                del self._successors[frozen_tail]
            # Remove frozen_tail as a predecessor of frozen_head
            del self._predecessors[frozen_head][frozen_tail]
            # If that head is no longer the head of any hyperedge, remove it
            # from the predecessors dictionary
            if self._predecessors[frozen_head] == {}:
                del self._predecessors[frozen_head]

        # Remove hyperedge's attributes dictionary
        del self._hyperedge_attributes[hyperedge_id]

    def remove_hyperedges(self, hyperedge_ids):
        """Removes a set of hyperedges and their attributes from
        the hypergraph.

        :param hyperedge_ids: iterable container of IDs of the hyperedges
                        to be removed.
        :raises: ValueError -- No such hyperedge exists.

        See also:
        remove_hyperedge

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> hyperedge_list = ((["A"], ["B", "C"]),
                                  (("A", "B"), ("C"), {'weight': 2}),
                                  (set(["B"]), set(["A", "C"])))
            >>> hyperedge_ids = H.add_hyperedges(hyperedge_list)
            >>> H.remove_hyperedges(hyperedge_ids)

        """
        for hyperedge_id in hyperedge_ids:
            self.remove_hyperedge(hyperedge_id)

    def has_undirected_hyperedge(self, nodes):
        """Given a set of nodes, returns whether there is a hyperedge in the
        hypergraph that is precisely composed of those nodes.

        :param nodes: iterable container of references to nodes in the
                    hyperedge being checked.
        :returns: bool -- true iff a hyperedge exists composed of the
                specified nodes.
        """
        frozen_nodes = frozenset(nodes)
        return frozen_nodes in self._node_set_to_hyperedge

    def has_directed_hyperedge(self, tail, head):
        """Given a tail and head set of nodes, returns whether there
        is a hyperedge in the hypergraph that connects the tail set
        to the head set.

        :param tail: iterable container of references to nodes in the
                    tail of the hyperedge being checked.
        :param head: iterable container of references to nodes in the
                    head of the hyperedge being checked.
        :returns: bool -- true iff a hyperedge exists connecting the
                specified tail set to the specified head set.
        """
        frozen_tail = frozenset(tail)
        frozen_head = frozenset(head)
        return frozen_tail in self._successors and \
            frozen_head in self._successors[frozen_tail]

    def has_hyperedge_id(self, hyperedge_id):
        """Determines if a hyperedge referenced by hyperedge_id
        exists in the hypergraph.

        :param hyperedge_id: ID of the hyperedge whose existence is
                            being checked.
        :returns: bool -- true iff a hyperedge exists that has id hyperedge_id.

        """
        return hyperedge_id in self._hyperedge_attributes

    def get_hyperedge_id_set(self):
        """Returns the set of IDs of hyperedges that are currently
        in the hypergraph.

        :returns: set -- all IDs of hyperedges currently in the hypergraph

        """
        return set(self._hyperedge_attributes.keys())

    def hyperedge_id_iterator(self):
        """Provides an iterator over the list of hyperedge IDs.

        """
        return iter(self._hyperedge_attributes)

    def get_undirected_hyperedge_id(self, nodes):
        """From a set of nodes, returns the ID of the undirected hyperedge
        that this set comprises.

        :param nodes: iterable container of references to nodes in the
                    the hyperedge to be added
        :returns: str -- ID of the hyperedge that has that the specified
                node set comprises.
        :raises: ValueError -- No such hyperedge exists.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> hyperedge_list = (["A", "B", "C"],
                                  ("A", "D"),
                                  set(["B", "D"]))
            >>> hyperedge_ids = H.add_undirected_hyperedges(hyperedge_list)
            >>> x = H.get_undirected_hyperedge_id(["A", "B", "C"])

        """
        frozen_nodes = frozenset(nodes)

        if not self.has_undirected_hyperedge(frozen_nodes):
            raise ValueError("No such hyperedge exists.")

        return self._node_set_to_hyperedge[frozen_nodes]

    def get_directed_hyperedge_id(self, tail, head):
        """From a tail and head set of nodes, returns the ID of the directed
        hyperedge that these sets comprise.

        :param tail: iterable container of references to nodes in the
                    tail of the hyperedge to be added
        :param head: iterable container of references to nodes in the
                    head of the hyperedge to be added
        :returns: str -- ID of the hyperedge that has that the specified
                tail and head sets comprise.
        :raises: ValueError -- No such hyperedge exists.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> hyperedge_list = (["A"], ["B", "C"]),
                                  (("A", "B"), ("C"), {weight: 2}),
                                  (set(["B"]), set(["A", "C"])))
            >>> hyperedge_ids = H.add_directed_hyperedges(hyperedge_list)
            >>> x = H.get_directed_hyperedge_id(["A"], ["B", "C"])

        """
        frozen_tail = frozenset(tail)
        frozen_head = frozenset(head)

        if not self.has_directed_hyperedge(frozen_tail, frozen_head):
            raise ValueError("No such hyperedge exists.")

        return self._successors[frozen_tail][frozen_head]

    def get_hyperedge_attribute(self, hyperedge_id, attribute_name):
        """Given a hyperedge ID and the name of an attribute, get a copy
        of that hyperedge's attribute.

        :param hyperedge_id: ID of the hyperedge to retrieve the attribute of.
        :param attribute_name: name of the attribute to retrieve.
        :returns: attribute value of the attribute_name key for the
                specified hyperedge.
        :raises: ValueError -- No such hyperedge exists.
        :raises: ValueError -- No such attribute exists.

        Examples:
        ::

            >>> H = MixedHypergraph()
            >>> hyperedge_list = (["A"], ["B", "C"]),
                                  (("A", "B"), ("C"), {weight: 2}),
                                  (set(["B"]), set(["A", "C"])))
            >>> hyperedge_ids = H.add_directed_hyperedges(hyperedge_list)
            >>> attribute = H.get_hyperedge_attribute(hyperedge_ids[0])

        """
        if not self.has_hyperedge_id(hyperedge_id):
            raise ValueError("No such hyperedge exists.")
        elif attribute_name not in self._hyperedge_attributes[hyperedge_id]:
            raise ValueError("No such attribute exists.")
        else:
            return copy.\
                copy(self._hyperedge_attributes[hyperedge_id][attribute_name])

    def get_hyperedge_attributes(self, hyperedge_id):
        """Given a hyperedge ID, get a dictionary of copies of that hyperedge's
        attributes.

        :param hyperedge_id: ID of the hyperedge to retrieve the attributes of.
        :returns: dict -- copy of each attribute of the specified hyperedge_id
                (except the private __frozen_tail and __frozen_head entries).
        :raises: ValueError -- No such hyperedge exists.

        """
        if not self.has_hyperedge_id(hyperedge_id):
            raise ValueError("No such hyperedge exists.")
        dict_to_copy = self._hyperedge_attributes[hyperedge_id].items()
        attributes = {}
        for attr_name, attr_value in dict_to_copy:
            if attr_name not in ("__frozen_nodes", "__frozen_tail", "__frozen_head"):
                attributes[attr_name] = copy.copy(attr_value)
        return attributes

    def get_directed_hyperedge_tail(self, hyperedge_id):
        """Given a directed hyperedge ID,
        get a copy of that hyperedge's tail.

        :param hyperedge_id: ID of the hyperedge to retrieve the tail from.
        :returns: a copy of the container of nodes that the user provided
                as the tail to the hyperedge referenced as hyperedge_id.

        """
        return self.get_hyperedge_attribute(hyperedge_id, "tail")

    def get_directed_hyperedge_head(self, hyperedge_id):
        """Given a directed hyperedge ID,
        get a copy of that hyperedge's head.

        :param hyperedge: ID of the hyperedge to retrieve the head from.
        :returns: a copy of the container of nodes that the user provided
            as the head to the hyperedge referenced as hyperedge_id.

        """
        return self.get_hyperedge_attribute(hyperedge_id, "head")

    def get_undirected_hyperedge_nodes(self, hyperedge_id):
        """Given an undirected hyperedge ID,
        get a copy of that hyperedge's nodes. This is equivalent to
        get_hyperedge_nodes, except with a runtime guarantee
        that the hyperedge_id is undirected.

        :param hyperedge_id: ID of the hyperedge to retrieve the nodes from.
        :returns: a copy of the container of nodes that the user provided
                for the hyperedge referenced as hyperedge_id.

        """
        return self.get_hyperedge_attribute(hyperedge_id, "nodes")

    def get_hyperedge_nodes(self, hyperedge_id):
        """Given a hyperedge ID, get a set of that hyperedge's nodes.
        If the hyperedge is undirected, it returns the underlying
        nodes. Otherwise, it returns the union of the head and tail nodes.

        :param hyperedge_id: ID of the hyperedge to retrieve the nodes from.
        :returns: a set of nodes or the tail union head, depending on the
                hyperedge, that the user provided for the hyperedge referenced
                as hyperedge_id.

        """
        if self.is_hyperedge_id_directed(hyperedge_id):
            head = self.get_directed_hyperedge_head(hyperedge_id)
            tail = self.get_directed_hyperedge_tail(hyperedge_id)
            return set(head).union(tail)
        else:
            return self.get_hyperedge_attribute(hyperedge_id, "nodes")

    def get_hyperedge_weight(self, hyperedge_id):
        """Given a hyperedge ID, get that hyperedge's weight.

        :param hyperedge: ID of the hyperedge to retrieve the weight from.
        :returns: a the weight of the hyperedge referenced as hyperedge_id.

        """
        return self.get_hyperedge_attribute(hyperedge_id, "weight")

    def get_forward_star(self, node):
        """Given a node, get a copy of that node's forward star related
        to directed hyperedges.

        :param node: node to retrieve the forward-star of.
        :returns: set -- set of hyperedge_ids for the hyperedges
                        in the node's forward star.
        :raises: ValueError -- No such node exists.

        """
        if node not in self._node_attributes:
            raise ValueError("No such node exists.")
        return self._forward_star[node].copy()

    def get_backward_star(self, node):
        """Given a node, get a copy of that node's backward star related
        to directed hyperedges.

        :param node: node to retrieve the backward-star of.
        :returns: set -- set of hyperedge_ids for the hyperedges
                in the node's backward star.
        :raises: ValueError -- No such node exists.

        """
        if node not in self._node_attributes:
            raise ValueError("No such node exists.")
        return self._backward_star[node].copy()

    def get_star(self, node):
        """Given a node, get a copy of that node's star, that is, the set of
        undirected hyperedges that the node belongs to.

        :param node: node to retrieve the star of.
        :returns: set -- set of hyperedge_ids for the undirected hyperedges
                        in the node's star.
        :raises: ValueError -- No such node exists.

        """
        if node not in self._node_attributes:
            raise ValueError("No such node exists.")
        return self._star[node].copy()

    def get_successors(self, tail):
        """Given a tail set of nodes, get a list of edges of which the node
        set is the tail of each edge.

        :param tail: set of nodes that correspond to the tails of some
                        (possibly empty) set of edges.
        :returns: set -- hyperedge_ids of the hyperedges that have tail
                in the tail.

        """
        frozen_tail = frozenset(tail)
        # If this node set isn't any tail in the hypergraph, then it has
        # no successors; thus, return an empty list
        if frozen_tail not in self._successors:
            return set()

        return set(self._successors[frozen_tail].values())

    def get_predecessors(self, head):
        """Given a head set of nodes, get a list of edges of which the node set
        is the head of each edge.

        :param head: set of nodes that correspond to the heads of some
                        (possibly empty) set of edges.
        :returns: set -- hyperedge_ids of the hyperedges that have head
                in the head.
        """
        frozen_head = frozenset(head)
        # If this node set isn't any head in the hypergraph, then it has
        # no predecessors; thus, return an empty list
        if frozen_head not in self._predecessors:
            return set()

        return set(self._predecessors[frozen_head].values())

    def copy(self):
        """Creates a new MixedHypergraph object with the same node and
        hyperedge structure.
        Copies of the nodes' and hyperedges' attributes are stored
        and used in the new hypergraph.

        :returns: MixedHypergraph -- a new hypergraph that is a copy of
                the current hypergraph

        """
        return self.__copy__()

    def __copy__(self):
        """Creates a new MixedHypergraph object with the same node and
        hyperedge structure.
        Copies of the nodes' and hyperedges' attributes are stored
        and used in the new hypergraph.

        :returns: MixedHypergraph -- a new hypergraph that is a copy of
                the current hypergraph

        """
        new_H = MixedHypergraph()

        # Loop over every node and its corresponding attribute dict
        # in the original hypergraph's _node_attributes dict
        for node, attr_dict in self._node_attributes.items():
            # Create a new dict for that node to store that node's attributes
            new_H._node_attributes[node] = {}
            # Loop over each attribute of that node in the original hypergraph
            # and, for each key, copy the corresponding value into the new
            # hypergraph's dictionary using the same key
            for attr_name, attr_value in attr_dict.items():
                new_H._node_attributes[node][attr_name] = \
                    copy.copy(attr_value)

        # Loop over every hyperedge_id and its corresponding attribute dict
        # in the original hypergraph's _hyperedge_attributes dict
        for hyperedge_id, attr_dict in self._hyperedge_attributes.items():
            # Create a new dict for that node to store that node's attributes
            new_H._hyperedge_attributes[hyperedge_id] = {}
            # Loop over each attribute of that hyperedge in the original
            # hypergraph and, for each key, copy the corresponding value
            # the new hypergraph's dictionary
            for attr_name, attr_value in attr_dict.items():
                new_H.\
                    _hyperedge_attributes[hyperedge_id][attr_name] = \
                    copy.copy(attr_value)

        # Copy all of the original hypergraph's stars 
        new_H._backward_star = self._backward_star.copy()
        for node in self._node_attributes.keys():
            new_H._backward_star[node] = \
                self._backward_star[node].copy()
            new_H._forward_star[node] = \
                self._forward_star[node].copy()
        new_H._star = self._star.copy()
        for node in self._node_attributes.keys():
            new_H._star[node] = self._star[node].copy()

        # Copy the original hypergraph's composed hyperedges
        for frozen_nodes, hyperedge_id in self._node_set_to_hyperedge.items():
            new_H._node_set_to_hyperedge[frozen_nodes] = \
                copy.copy(hyperedge_id)
        # Copy the original hypergraph's successors
        for frozen_tail, successor_dict in self._successors.items():
            new_H._successors[frozen_tail] = successor_dict.copy()
        # Copy the original hypergraph's predecessors
        for frozen_head, predecessor_dict in self._predecessors.items():
            new_H._predecessors[frozen_head] = predecessor_dict.copy()

        # Start assigning edge labels at the same
        new_H._current_hyperedge_id = self._current_hyperedge_id

        return new_H

    def get_symmetric_image(self):
        """Creates a new MixedHypergraph object that is the symmetric
        image of this hypergraph (i.e., identical hypergraph with all
        edge directions reversed).
        Copies of each of the nodes' and hyperedges' attributes are stored
        and used in the new hypergraph.
        This does not affect undirected hyperedges.

        :returns: MixedHypergraph -- a new hypergraph that is the symmetric
                image of the current hypergraph.

        """
        new_H = self.copy()

        # No change to _node_attributes necessary, as nodes remain the same

        # Reverse the tail and head (and __frozen_tail and __frozen_head) for
        # every hyperedge
        for hyperedge_id in self.get_hyperedge_id_set():
            if self.is_hyperedge_id_directed(hyperedge_id):
                attr_dict = new_H._hyperedge_attributes[hyperedge_id]
                attr_dict["tail"], attr_dict["head"] = \
                    attr_dict["head"], attr_dict["tail"]
                attr_dict["__frozen_tail"], attr_dict["__frozen_head"] = \
                    attr_dict["__frozen_head"], attr_dict["__frozen_tail"]

        # Reverse the definition of forward star and backward star
        new_H._forward_star, new_H._backward_star = \
            new_H._backward_star, new_H._forward_star

        # Reverse the definition of successor and predecessor
        new_H._successors, new_H._predecessors = \
            new_H._predecessors, new_H._successors

        return new_H

    def get_induced_subhypergraph(self, nodes):
        """Gives a new hypergraph that is the subhypergraph of the current
        hypergraph induced by the provided set of nodes. That is, the induced
        subhypergraph's node set corresponds precisely to the nodes provided,
        and the coressponding hyperedges in the subhypergraph are only those
        from the original graph consist of tail and head sets that are subsets
        of the provided nodes.

        :param nodes: the set of nodes to find the induced subhypergraph of.
        :returns: MixedHypergraph -- the subhypergraph induced on the
                provided nodes.

        """
        sub_H = self.copy()
        sub_H.remove_nodes(sub_H.get_node_set() - set(nodes))
        return sub_H

    # TODO: make reading more extensible (attributes, variable ordering, etc.)
    def read(self, file_name, delim=',', sep='\t'):
        """Read a hypergraph from a file, where nodes are
        represented as strings.
        Each column is separated by "sep", and the individual
        tail nodes and head nodes are delimited by "delim".
        The header line is currently ignored, but columns should be of
        the format:
        tailnode1[delim]..tailnodeM[sep]headnode1[delim]..headnodeN[sep]direction[sep]weight

        Direction is required and must be (case-insensitive) U for undirected
        hyperedge or D for directed hyperedge. In the case of an undirected hyperedge,
        the 'tail' and 'head' nodes are added to the same nodes set to be added.
        
        As a concrete example, an arbitrary line with delim=',' and
        sep='    ' (4 spaces) may look like:
        ::

            x1,x2    x3,x4,x5    D    12

        which defines a directed hyperedge of weight 12 from a tail set containing
        nodes "x1" and "x2" to a head set containing nodes "x3", "x4", and "x5."

        An undirected hyperedge can also leave the head or tail set empty. As an example,
        the following undirected hyperedge representations are equivalent:
        ::

            x1    x2,x3    U
            x1,x2    x2,x3    U
            x1,x2,x3        U
                x1,x2,x3    U

        To read in a specific direction instead, use the other hypergraph type first,
        then convert it into a mixed hypergraph.
        """
        in_file = open(file_name, 'r')

        # Skip the header line
        in_file.readline()

        line_number = 2
        for line in in_file.readlines():
            line = line.rstrip()
            # Skip empty lines
            if not line:
                continue

            words = line.split(sep)
            if not (3 <= len(words) <= 4):
                raise \
                    IOError("Line {} ".format(line_number) +
                            "contains {} ".format(len(words)) +
                            "columns -- must contain only 3 or 4.")

            tail = set([word for word in words[0].split(delim)  
                if not word.isspace() and not word == ''])
            head = set([word for word in words[1].split(delim)  
                if not word.isspace() and not word == ''])

            direction = words[2]
            if len(words) == 4:
                weight = float(words[3].split(delim)[0])
            else:
                weight = 1
            if direction.upper() == 'D':
                self.add_directed_hyperedge(tail, head, weight=weight)
            elif direction.upper() == 'U':
                self.add_undirected_hyperedge(tail | head, weight=weight)
            else:
                raise IOError("Line {} ".format(line_number) +
                              "has direction {}, ".format(direction) +
                              "but it must be either U or D, " +
                              "case insensitive.")

            line_number += 1

        in_file.close()

    # TODO: make writing more extensible (attributes, variable ordering, etc.)
    def write(self, file_name, delim=',', sep='\t'):
        """Write a mixed hypergraph to a file, where nodes are
        represented as strings.
        Each column is separated by "sep", and the individual
        tail nodes and head nodes are delimited by "delim".
        The header line is currently ignored, but columns should be of
        the format:
        tailnode1[delim]..tailnodeM[sep]headnode1[delim]..headnodeN[sep]direction[sep]weight

        As a concrete example, an arbitrary line with delim=',' and
        sep='    ' (4 spaces) may look like:
        ::

            x1,x2    x3,x4,x5    D    12

        which defines a directed hyperedge of weight 12 from a tail set containing
        nodes "x1" and "x2" to a head set containing nodes "x3", "x4", and "x5",

        For undirected hyperedges, `write` prefers to add the entire `nodes` set to the
        tail column instead. For example:
        ::

            x1,x2,x3        U    15

        """
        out_file = open(file_name, 'w')

        # write first header line
        out_file.write("tail" + sep + "head" + sep + "direction" + sep + "weight\n")

        for hyperedge_id in self.get_hyperedge_id_set():
            line = ""
            if self.is_hyperedge_id_undirected(hyperedge_id):
                # Write each node to the line, separated by delim
                for node in self.get_undirected_hyperedge_nodes(hyperedge_id):
                    line += node + delim
                # Remove last (extra) delim
                line = line[:-1]
                # Add extra sep to account for no head nodes
                line += sep
            else:
                # Write each tail node to the line, separated by delim
                for tail_node in self.get_directed_hyperedge_tail(hyperedge_id):
                    line += tail_node + delim
                # Remove last (extra) delim
                line = line[:-1]

                line += sep

                # Write each head node to the line, separated by delim
                for head_node in self.get_directed_hyperedge_head(hyperedge_id):
                    line += head_node + delim
                # Remove last (extra) delim
                line = line[:-1]

            # Add sep between columns
            line += sep

            # Add the corresponding direction letter
            if self.is_hyperedge_id_undirected(hyperedge_id):
                line += "U"
            else:
                line += "D"

            # Write the weight to the line and end the line
            line += sep + str(self.get_hyperedge_weight(hyperedge_id)) + "\n"

            out_file.write(line)

        out_file.close()

    def underlying_undirected_hypergraph(self):
        """Constructs an underlying undirected hypergraph,
        which contains all the nodes of this mixed hypergraph
        plus all of its undirected hyperedges.

        :returns: a backing UndirectedHypergraph
        """
        undirected_hypergraph = UndirectedHypergraph()

        for node in self.get_node_set():
            undirected_hypergraph.add_node(node)
        
        for hyperedge_id in self.get_hyperedge_id_set():
            if self.is_hyperedge_id_undirected(hyperedge_id):
                nodes = self.get_undirected_hyperedge_nodes(hyperedge_id)
                attrs = self.get_hyperedge_attributes(hyperedge_id)
                undirected_hypergraph.add_hyperedge(nodes, attrs)
        
        return undirected_hypergraph
    
    def underlying_directed_hypergraph(self):
        """Constructs an underlying directed hypergraph,
        which contains all the nodes of this mixed hypergraph
        plus all of its directed hyperedges.

        :returns: a backing DirectedHypergraph
        """
        directed_hypergraph = DirectedHypergraph()

        for node in self.get_node_set():
            directed_hypergraph.add_node(node)
        
        for hyperedge_id in self.get_hyperedge_id_set():
            if self.is_hyperedge_id_directed(hyperedge_id):
                tail = self.get_directed_hyperedge_tail(hyperedge_id)
                head = self.get_directed_hyperedge_head(hyperedge_id)
                attrs = self.get_hyperedge_attributes(hyperedge_id)
                directed_hypergraph.add_hyperedge(tail, head, attrs)
        
        return directed_hypergraph
    
    def extend_directed_hypergraph(self, directed_hypergraph):
        """Appends all of the nodes and directed hyperedges
        of the input directed hypergraph to this mixed hypergraph.
        """
        for node in directed_hypergraph.get_node_set():
            self.add_node(node)
        
        for hyperedge_id in directed_hypergraph.get_hyperedge_id_set():
            tail = directed_hypergraph.get_hyperedge_tail(hyperedge_id)
            head = directed_hypergraph.get_hyperedge_head(hyperedge_id)
            attrs = directed_hypergraph.get_hyperedge_attributes(hyperedge_id)
            self.add_directed_hyperedge(tail, head, attrs)
    
    def extend_undirected_hypergraph(self, undirected_hypergraph):
        """Appends all of the nodes and undirected hyperedges
        of the input undirected hypergraph to this mixed hypergraph.
        """
        for node in undirected_hypergraph.get_node_set():
            self.add_node(node)
        
        for hyperedge_id in undirected_hypergraph.get_hyperedge_id_set():
            nodes = undirected_hypergraph.get_hyperedge_nodes(hyperedge_id)
            attrs = undirected_hypergraph.get_hyperedge_attributes(hyperedge_id)
            self.add_undirected_hyperedge(nodes, attrs)
