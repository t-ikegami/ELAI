//
// ELAI User's Custom Class-Boilerplate ( EUCC-B )
//

class Element
{
  // DATA IMPLEMENTATION

public:
  typedef /* USER DEFINED */ id_type;

  Element( const Element& x )
  {}
  Element( const Element& x, int c ) // FOR COLORING
  {}
  Element( const id_type& )
  {}
  ~Element() {}

  id_type operator()() const
  {}

  int color() const
  {}

  bool operator<( const Element& x ) const
  {}
  bool operator==( const Element& x ) const
  {}
};

class Neighbour
{
  // DATA IMPLEMENTATION

public:
  typedef /* IMPLEMENTED CONTAINER */::iterator iterator;
  typedef /* SAME ABOVE */::const_iterator const_iterator;

  Neighbour( const Element& x )
  {}
  Neighbour( const Neighbour& u )
  {}
  ~Neighbour() {}

  Neighbour& operator=( const Neighbour& src )
  {}

  const Element& element() const
  {}

  int size() const
  {}

  Neighbour& join( const Element& x )
  {}

  Neighbour& erase( const Element& x )
  {}
  Neighbour& erase( const_iterator& it )
  {}

  bool operator()( const Element& x ) const
  {}

  iterator begin()
  {}
  const_iterator begin() const
  {}
  const_iterator end() const
  {}
};

