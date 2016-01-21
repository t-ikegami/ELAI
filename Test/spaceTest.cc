#include <iostream>
#include "space.hpp"

using namespace std;

class Element
{
public:
  enum type { A, B, C };

private:
  const int id;
  const type t;

public:
  Element( int id0, type t0 ) : id( id0 ), t( t0 ) {}
  Element( const Element& src ) : id( src.id ), t( src.t ) {}
  ~Element() {}

  bool operator<( const Element& rhs ) const
  {
    if ( id == rhs.id ) return t < rhs.t;
    return id < rhs.id;
  }
  bool operator==( const Element& rhs ) const
  {
    return id == rhs.id && t == rhs.t;
  }
};

typedef elai::space< Element > Space;

int main()
{
  {
    Space s1, s2, s3;

    s1.join( Element( 0, Element::A ) );
    s1.join( Element( 0, Element::B ) );
    s1.join( Element( 0, Element::C ) );
    s1.join( Element( 1, Element::A ) );
    s1.join( Element( 2, Element::A ) );
    s1.join( Element( 2, Element::B ) );
    s1.join( Element( 2, Element::C ) );

    s2.join( Element( 0, Element::A ) );
    s2.join( Element( 1, Element::A ) );
    s2.join( Element( 2, Element::A ) );

    s3.join( Element( 0, Element::A ) );
    s3.join( Element( 1, Element::A ) );
    s3.join( Element( 2, Element::A ) );

    s3 = s1 & s2;
    assert( s3 == s2 );
  }

  {
    Space s1, s2, s3;

    s1.join( Element( 0, Element::A ) );
    s1.join( Element( 1, Element::A ) );
    s1.join( Element( 2, Element::A ) );

    s2.join( Element( 0, Element::B ) );
    s2.join( Element( 2, Element::B ) );

    s3.join( Element( 0, Element::A ) );
    s3.join( Element( 0, Element::B ) );
    s3.join( Element( 1, Element::A ) );
    s3.join( Element( 2, Element::A ) );
    s3.join( Element( 2, Element::B ) );

    assert( s3 == ( s1 | s2 ) );
  }

  {
    Space s1, s2, s3;

    s1.join( Element( 0, Element::A ) );
    s1.join( Element( 0, Element::B ) );
    s1.join( Element( 0, Element::C ) );
    s1.join( Element( 1, Element::A ) );
    s1.join( Element( 2, Element::A ) );
    s1.join( Element( 2, Element::B ) );
    s1.join( Element( 2, Element::C ) );

    s2.join( Element( 0, Element::A ) );
    s2.join( Element( 1, Element::A ) );
    s2.join( Element( 2, Element::A ) );

    s3.join( Element( 0, Element::B ) );
    s3.join( Element( 0, Element::C ) );
    s3.join( Element( 2, Element::B ) );
    s3.join( Element( 2, Element::C ) );

    s1 /= s2;
    assert( s3 == s1 );
  }
}
