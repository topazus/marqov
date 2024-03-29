/***************************************************************************
 *   Copyright (C) 2005-2022 by Florian Goth   *
 *   fgoth@physik.uni-wuerzburg.de   *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#ifndef REGISTRY_H
#define REGISTRY_H

/** @file Registry include header.
 *
 * This includes the entry point class, the exception types and some auxiliary
 * functions.
 * @author Florian Goth.
 *
 */
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>

/** The base exception class thrown by the registry.
 *
 * It is derived from the STL logic_error exception.
 */
class Registry_Exception : public std::logic_error
{
public:
    /** Construct exception with error message.
     * 
     * @param err_msg The error message.
     */
    explicit Registry_Exception(const std::string& err_msg) throw() : std::logic_error(err_msg) {}
private:
};

/** The exception when the desired key is not found.
 */
class Registry_Key_not_found_Exception : public Registry_Exception
{
public:
    /** The Exception when a requested key is not found.
     * 
     * @param err_msg The error message of which key is not found.
     */
    explicit Registry_Key_not_found_Exception(const std::string& err_msg) throw() : Registry_Exception(err_msg) {}
private:
};

/** The exception when the requested data in a block is not present.
 */
class Registry_Block_Data_not_found_Exception : public Registry_Key_not_found_Exception
{
public:
    const std::string block_key; ///< to denote the block where no data was found
    /** Construct an exception when the data was not found.
     * 
     * @param key_value which key was not found.
     */
    explicit Registry_Block_Data_not_found_Exception(const std::string& key_value) throw() : Registry_Key_not_found_Exception(std::string("Block Key not found: ") + key_value), block_key(key_value) {}
    ~Registry_Block_Data_not_found_Exception() throw() {}
private:
};

/** The exception when the requested block is not present.
 */
class Registry_Block_not_found_Exception : public Registry_Key_not_found_Exception
{
public:
    const std::string block;///< which block has not been found.
    /** Construct an exception where a block was not found.
     * 
     * @param key_value the block which has not been found.
     */
    explicit Registry_Block_not_found_Exception(const std::string& key_value) throw() : Registry_Key_not_found_Exception(std::string("Block not found: ") + key_value), block(key_value) {}
    ~Registry_Block_not_found_Exception() throw() {}
private:
};

/** The exception when the requested config file could not be found.
 */
class Registry_cfgfile_not_found_Exception : public Registry_Key_not_found_Exception
{
public:
    const std::string cfgfile; ///< which config file was not found.
    /** Construct the respective exception when a file was not found.
     * 
     * @param key_value Which config file could not be found.
     */
    explicit Registry_cfgfile_not_found_Exception(const std::string& key_value) throw(): Registry_Key_not_found_Exception(std::string("cfgfile not found: ") + key_value), cfgfile(key_value) {}
    ~Registry_cfgfile_not_found_Exception() throw() {}
private:
};

/**
 * This class contains the contents of a [BLOCK] structure in a config file.
 * 
 * In Config Files those [BLOCK] thingies.
 */
class Block
{
    std::map<std::string , std::string > Block_Data;///< A map between the keys and their values.
    std::string BlockName; ///< The name of the block
    std::vector<std::string> Keys; ///< all keys that have been found.
    std::map<std::string , std::string >::size_type NrOfKeys; ///< the number of key-value pairs that we have.
    /** Add a key-value pair.
     * 
     * @param arg Add this key-value pair.
     */
    void push_back_Key(std::string arg);
public:
    /** Create an empty block
     */
    Block() {}
    /** Create a block from a vector of strings.
     * 
     * @param blockname The name of the block.
     * @param block a vector of strings.
     */
    Block(std::string blockname , const std::vector<std::string>& block);
    /** Get the number of keys that we store.
     * 
     * @returns the number of keys that we store.
     */
    inline unsigned int size(void) const;
    /** Get the list of keys that we store
     * 
     * @returns The list of keys that we store.
     */
    const std::vector<std::string> GetKeys() const
    {
        return Keys;
    }
    /** Get access to the contents of a block.
     * 
     * This distinguishes itself from find() by doing no checks.
     * @param key the key for which to look.
     * @returns the value of the key.
     */
    std::string& operator[](const std::string& key)
    {
        return Block_Data[key];
    }
    /** Get access to the contents of a block.
     * 
     * non-const version.
     * In contrast to operator[] this perform checking.
     * @param key the key for which to look.
     * @returns the value of the key.
     * @throws Registry_Block_Data_not_found_Exception If the key was not found
     */
    std::string& find(const std::string& key)
    {
        std::map<std::string , std::string >::iterator it = Block_Data.find(key);
        if (it == Block_Data.end()) throw(Registry_Block_Data_not_found_Exception(key));
        return it->second;
    }
    /** Get access to the contents of a block.
     * 
     * constified version.
     * In contrast to operator[] this perform checking.
     * @param key the key for which to look.
     * @returns the value of the key.
     * @throws Registry_Block_Data_not_found_Exception If the key was not found
     */
    const std::string& find(const std::string& key) const
    {
        std::map<std::string , std::string >::const_iterator it = Block_Data.find(key);
        if (it == Block_Data.end()) throw(Registry_Block_Data_not_found_Exception(key));
        return it->second;
    }
    /** Clean up a block.
     */
    ~Block()
    {}
};

/**
 * This class contains the contents of a single config file.
 */
class cfile
{
    std::string filename; ///< the filename
    std::map<std::string , Block> cfgfile; ///< The map between files and their internal blocks.
public:
    /** Default constructor.
     */
    cfile() {}
    /** Get access to the contents of a file.
     * 
     * This distinguishes itself from find() by doing no checks.
     * @param a the block for which to look.
     * @returns A block that can be queried further.
     */
    Block& operator[](std::string a)
    {
        return cfgfile[a];
    }
    /** Get access to the contents of a file.
     * 
     * non-const version.
     * In contrast to operator[] this perform checking.
     * @param key the file for which to look.
     * @returns the blocks of the file.
     */
    Block& find(const std::string& key)
    {
        std::map<std::string , Block>::iterator it = cfgfile.find(key);
        if ( it == cfgfile.end()) throw(Registry_Block_not_found_Exception(key));
        return it->second;
    }
    /** Get access to the contents of a file.
     * 
     * const version.
     * In contrast to operator[] this perform checking.
     * @param key the file for which to look.
     * @returns the blocks of the file.
     */
    const Block& find(const std::string& key) const
    {
        std::map<std::string , Block>::const_iterator it = cfgfile.find(key);
        if ( it == cfgfile.end()) throw(Registry_Block_not_found_Exception(key));
        return it->second;
    }
    /**
     * Constructor that initializes this with the contents of a config-file.
     * 
     * @param file The file that we should parse.
     */
    cfile(std::string& file);
};

/**
 * This holds together all the contents of the configuration directory and
 * provides access via the Get() template.
 */
class RegistryDB
{
    std::map<std::string , cfile> Reg; ///< The map containing all config files.
private:
public:
    /** Initialize registry.
     * 
     * This initializes the registry via a separate function call.
     * @param cfgDir the directory that contains all the files the registry should contain
     * @param pat a suffix to select only certain files, e.g. : .ini
     */
    int init(const std::string& cfgDir, const std::string pat = "");
    /** construct the registry.
     * 
     * This is the constructor for the registry.
     * @param arg the directory that contains all the files the registry should contain
     * @param pat a suffix to select only certain files, e.g. : .ini
     */
    RegistryDB(const std::string& arg, const std::string pat = "");
    /** Empty default constructor.
     */
    RegistryDB()
    {}
    /** Tidy up everything.
     */
    ~RegistryDB()
    {}
    /** Get a particular block in a file of the Registry.
     * 
     * @param file the file in which to look.
     * @param bloc the bloc we require.
     * @returns the requested block.
     */
    Block GetBlock(std::string file , std::string bloc)
    {
        return Reg[file][bloc];
    }
    /** Function to get a value from the registry.
     *
     * The template parameter determines to which type to convert the key.
     * @tparam T to which type do we convert.
     * 
     * @param file the file in which to look
     * @param block under which block is the value
     * @param key for which key to look
     * @return the requested key
     */
    template < typename T >
    inline T Get(const std::string& file, const std::string& block, const std::string& key) const;
    /** A function for setting values in the registry.
     */
    template <typename T>
    inline T set(const std::string& file, const std::string& block, const std::string& key, T val);
};

/** The helper template for performing string -> type conversions.
 * 
 * The basic template for doing the conversion between strings and the requested type.
 * We use the C++ stringstreams thus we benefit from all overloads that are already provided by C++.
 * @tparam A to which type do we want to convert.
 */
template < typename A >
struct GetTrait
{
    /** Helper function to convert from a string to the requested type.
     * @param arg a string that should represent something.
     * @return Hopefully, the successfully converted object.
     */
    static inline A Convert(std::string arg)
    {
        A ret;
        std::stringstream sstream(arg);
        sstream >> ret;
        return ret;
    }
};

/** The helper template for performing string -> bool conversion
 * 
 * Specialization for boolean(true, false) like strings
 */
template <>
struct GetTrait<bool>
{
    /** Helper function to convert from a string to a boolean value
     * 
     * Every occurence of uppercase/lowercase mixing of TRUE is interpreted as true,
     * everything else is false.
     * @return boolean true, if the string was [Tt][Rr][Uu][Ee]
     */
    static bool Convert(std::string arg)
    {
        std::transform ( arg.begin(), arg.end(), arg.begin(), ::toupper );
        if ( arg == "TRUE")
            return true;
        return false;
    }
};

/** A helper template to set values in the registry.
 * 
 * @tparam A the type of the value.
 */
template < typename A >
struct SetTrait
{
    /** Convert a type to its textual representation.
     *
     * @param arg The value we want to write.
     * @return a textual representation of arg.
     */
    static std::string convert(A arg)
    {
        std::stringstream ss;
        ss << arg << '\0';
        std::string os(ss.str());
        return os;
    }
};

template < typename T >
inline T RegistryDB::Get(const std::string& file, const std::string& block, const std::string& Key) const
{
    typedef GetTrait<T> Trait;
    std::map<std::string , cfile>::const_iterator tempkey = Reg.find(file);
    if (tempkey == Reg.end()) throw(Registry_cfgfile_not_found_Exception(std::string("File Key not found :") + file));
    std::string temp(tempkey->second.find(block).find(Key));
    return Trait::Convert(temp);
}

template <typename T>
inline T RegistryDB::set(const std::string& file, const std::string& block, const std::string& key, T val)
{
    typedef SetTrait<T> Trait;
    Reg[file][block][key] = Trait::convert(val);
    return val;
}

inline unsigned int Block::size(void) const
{
    return Block_Data.size();
}

//Two examples on how to extend the parsing capabilities of the registry.

/** A helper trait for reading vectors of values with a predefined separator.
 * 
 * The predefined separator is currently hard-coded to ";" or ",".
 * @tparam T the type of the elements in the vector.
 */
template <typename T>
struct GetTrait<std::vector<T> >
{
    /** Implementation function for the conversion.
     * 
     * @param arg the string that we intend to break up.
     * @return the converted values stored in a C++ std::vector .
     */
    static std::vector<T> Convert(std::string& arg)
    {
        const std::string delim("; ,");
        std::vector<T> retval;
        std::size_t pos = 0;
        std::size_t posold = 0;
        do
        {
            retval.push_back(GetTrait<T>::Convert(arg.substr(posold, ( pos = arg.find_first_of(delim, posold) ) - posold) ));
        } while ((posold = arg.find_first_not_of(delim, pos) ) != std::string::npos);
        return retval;
    }
};

/** A helper trait for reading vectors of strings with a predefined separator.
 * 
 * The predefined separator is currently hard-coded to ";" or ",".
 * @tparam T the type of the elements in the vector.
 */
template <>
struct GetTrait<std::vector<std::string> >
{
    /** Implementation function for the conversion.
     * 
     * @param arg the string that we intend to break up.
     * @return the strings stored in a C++ std::vector .
     */
    static std::vector<std::string> Convert(std::string arg)
    {
        const std::string delim("; ,");
        std::vector<std::string> retval;
        std::size_t pos = 0;
        std::size_t posold = 0;
        do
        {
            retval.push_back(arg.substr(posold, ( pos = arg.find_first_of(delim, posold) ) - posold) );
        } while ((posold = arg.find_first_not_of(delim, pos) ) != std::string::npos);
        return retval;
    }
};

#endif
