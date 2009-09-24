/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines the Property and PropertySet classes
/// \todo Improve documentation


#ifndef __defined_libdai_properties_h
#define __defined_libdai_properties_h


#include <iostream>
#include <sstream>
#include <boost/any.hpp>
#include <map>
#include <vector>
#include <typeinfo>
#include <dai/exceptions.h>
#include <dai/util.h>
#include <boost/lexical_cast.hpp>


namespace dai {


/// Type of the key of a Property
typedef std::string PropertyKey;

/// Type of the value of a Property
typedef boost::any  PropertyValue;

/// A Property is a pair of a key and a corresponding value
typedef std::pair<PropertyKey, PropertyValue> Property;


/// Writes a Property object to an output stream
std::ostream& operator<< (std::ostream & os, const Property & p);


/// Represents a set of properties, mapping keys (of type PropertyKey) to values (of type PropertyValue)
class PropertySet : private std::map<PropertyKey, PropertyValue> {
    public:
        /// Default constructor
        PropertySet() {}

        /// Construct PropertySet from a string
        PropertySet( const std::string &s ) {
            std::stringstream ss;
            ss << s;
            ss >> *this;
        }

        /// Gets a property
        const PropertyValue & Get(const PropertyKey &key) const {
            PropertySet::const_iterator x = find(key);
            if( x == this->end() )
                DAI_THROWE(OBJECT_NOT_FOUND,"PropertySet::Get cannot find property '" + key + "'");
            return x->second;
        }

        /// Sets a property
        PropertySet & Set(const PropertyKey &key, const PropertyValue &val) { this->operator[](key) = val; return *this; }

        /// Set properties according to those in newProps, overriding properties that already exist with new values
        PropertySet & Set( const PropertySet &newProps ) {
            const std::map<PropertyKey, PropertyValue> *m = &newProps;
            foreach(value_type i, *m)
                Set( i.first, i.second );
            return *this;
        }

        /// Gets a property, casted as ValueType
        template<typename ValueType>
        ValueType GetAs(const PropertyKey &key) const {
            try {
                return boost::any_cast<ValueType>(Get(key));
            } catch( const boost::bad_any_cast & ) {
                DAI_THROWE(IMPOSSIBLE_TYPECAST,"Cannot cast value of property '" + key + "' to desired type.");
                return ValueType();
            }
        }

        /// Converts a property from string to ValueType (if necessary)
        template<typename ValueType>
        void ConvertTo(const PropertyKey &key) { 
            PropertyValue val = Get(key);
            if( val.type() != typeid(ValueType) ) {
                DAI_ASSERT( val.type() == typeid(std::string) );
                try {
                    Set(key, boost::lexical_cast<ValueType>(GetAs<std::string>(key)));
                } catch(boost::bad_lexical_cast &) {
                    DAI_THROWE(IMPOSSIBLE_TYPECAST,"Cannot cast value of property '" + key + "' from string to desired type.");
                }
            }
        }

        /// Converts a property from string to ValueType (if necessary)
        template<typename ValueType>
        ValueType getStringAs(const PropertyKey &key) const { 
            PropertyValue val = Get(key);
            if( val.type() == typeid(ValueType) ) {
                return boost::any_cast<ValueType>(val);
            } else if( val.type() == typeid(std::string) ) {
                try {
                    return boost::lexical_cast<ValueType>(GetAs<std::string>(key));
                } catch(boost::bad_lexical_cast &) {
                    DAI_THROWE(IMPOSSIBLE_TYPECAST,"Cannot cast value of property '" + key + "' from string to desired type.");
                }
            } else
                DAI_THROWE(IMPOSSIBLE_TYPECAST,"Cannot cast value of property '" + key + "' from string to desired type.");
            return ValueType();
        }

        /// Converts a property from ValueType to string (if necessary)
        template<typename ValueType>
        PropertySet & setAsString(const PropertyKey &key, ValueType &val) {
            try {
                return Set( key, boost::lexical_cast<std::string>(val) );
            } catch( boost::bad_lexical_cast & ) {
                DAI_THROWE(IMPOSSIBLE_TYPECAST,"Cannot cast value of property '" + key + "' to string.");
            }
        }

        /// Shorthand for (temporarily) adding properties, e.g. PropertySet p()("method","BP")("verbose",1)("tol",1e-9)
        PropertySet operator()(const PropertyKey &key, const PropertyValue &val) const { PropertySet copy = *this; return copy.Set(key,val); }

        /// Check if a property with the given key exists
        bool hasKey(const PropertyKey &key) const { PropertySet::const_iterator x = find(key); return (x != this->end()); }

        /// Returns a set containing all keys
        std::set<PropertyKey> allKeys() const {
            std::set<PropertyKey> res;
            const_iterator i;
            for( i = begin(); i != end(); i++ )
                res.insert( i->first );
            return res;
        }

        /// Returns a vector containing all keys
        std::vector<PropertyKey> keys() const {
            std::vector<PropertyKey> result;
            result.reserve( size() );
            for( PropertySet::const_iterator i = begin(); i != end(); ++i )
                result.push_back( i->first );
            return result;
        }

        /// Writes a PropertySet object to an output stream
        friend std::ostream& operator<< (std::ostream & os, const PropertySet & ps);

        /// Reads a PropertySet object from an input stream, storing values as strings
        friend std::istream& operator>> (std::istream& is, PropertySet & ps);
};


} // end of namespace dai


#endif
