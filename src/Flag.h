/** \file Flagh.h */

#pragma once

namespace AMDiS
{
  /** \ingroup Common
   * \brief
   * The Flag class encapsulates flags which represents simple information.
   * Used e.g. while mesh traversal to specify, which elements should be
   * visited.
   */
  class Flag
  {
  public:
    /// Constructs a unset Flag
    Flag() : flags(0) {}

    /// Constructs a Flag initialized by f
    Flag(const unsigned long f) : flags(f) {}

    /// Copy constructor
    Flag(Flag const& f) : flags(f.flags) {}

    /// Compares two Flags
    bool operator==(Flag const& f) const
    {
      return (flags == f.flags);
    }

    /// Compares two Flags
    bool operator!=(Flag const& f) const
    {
      return !(f == *this);
    }

    /// Assignment operator
    Flag& operator=(Flag const& f)
    {
      if (this != &f)
        flags = f.flags;
      return *this;
    }

    /// Typecast
    operator bool() const
    {
      return isAnySet();
    }

    /// Set \ref flags
    void setFlags(const unsigned long f)
    {
      flags = f;
    }

    /// Set \ref flags
    void setFlags(Flag const& f)
    {
      flags = f.flags;
    }

    /// Sets \ref flags to \ref flags | f
    void setFlag(const unsigned long f)
    {
      flags |= f;
    }

    /// Sets \ref flags to \ref flags | f.flags
    void setFlag(Flag const& f)
    {
      flags |= f.flags;
    }

    /// Sets \ref flags to \ref flags & ~f
    void unsetFlag(const unsigned long f)
    {
      flags &= ~f;
    }

    /// Sets \ref flags to \ref flags & ~f.flags
    void unsetFlag(Flag const& f)
    {
      flags &= ~f.flags;
    }

    unsigned long getFlags() const
    {
      return flags;
    }

    /// Returns \ref flags | f.flags
    Flag operator+(Flag const& f) const
    {
      Flag r(flags);
      r.setFlag(f);
      return r;
    }

    /// Returns \ref flags & ~f.flags
    Flag operator-(Flag const& f) const
    {
      Flag r(flags);
      r.unsetFlag(f);
      return r;
    }

    /// Returns \ref flags | f.flags
    Flag operator|(Flag const& f) const
    {
      Flag r(flags);
      r.setFlag(f);
      return r;
    }

    /// Returns \ref flags & f.flags
    Flag operator&(Flag const& f) const
    {
      Flag r(flags);
      r.flags &= f.flags;
      return r;
    }

    /// Sets \ref flags to \ref flags &= f.flags
    Flag operator&=(Flag const& f)
    {
      flags &= f.flags;
      return *this;
    }

    /// Returns \ref flags ^ f.flags
    Flag operator^(Flag const& f) const
    {
      Flag r(flags);
      r.flags ^= f.flags;
      return r;
    }

    /// Sets \ref flags to \ref flags & f.flags
    Flag& operator|=(Flag const& f)
    {
      if (this != &f)
        flags |= f.flags;
      return *this;
    }

    /// Returns ~\ref flags
    Flag operator~() const
    {
      Flag r;
      r.flags = ~flags;
      return r;
    }

    /// Checks whether all set bits of f.flags are set in \ref flags too.
    bool isSet(Flag const& f) const
    {
      return ((flags&f.flags) == f.flags);
    }

    /// Returns !\ref isSet(f)
    bool isUnset(Flag const& f) const
    {
      return !isSet(f);
    }

    /// Returns true if \ref flags != 0
    bool isAnySet() const
    {
      return (flags != 0);
    }

  protected:
    /// Internal flag representation
    unsigned long  flags;
  };

} // end namespace AMDiS
