//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_BASE_MAT4_HH
#define VC_BASE_MAT4_HH

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <map>

#include "../system/platform.hh" // alloca

#include "../base/endian.hh"
#include "../base/exception.hh"
#include "../base/printf.hh"
#include "../base/primtype.hh"

//=============================================================================

# if defined(_MSC_VER)
#  pragma warning (push)
#  pragma warning (disable: 4290)
# endif

namespace VC {
namespace base {

/** Read/write MATLAB level 4 MAT file format \ingroup vc_mat4

    This `namespace` is aliased as `VC::math::mat4`.
 */
namespace mat4io {

/** \defgroup vc_mat4 Read/write MATLAB level 4 MAT file format
    \ingroup vc_base

    MATLAB's level 4 format is outdated: you have to force MATLAB
    writing MAT4 files, e.g., "save 'myfile.mat' -V4".

    Why is this cool, though?
    - The level 4 file format is simple.
    - It's good enough to exchange data easily with MATLAB/Octave/...
    - It's fine for streaming.

    The last point means: MAT4 encodes each matrix separately. There is
    no directory structure, no hierarchy (and you can store essentially
    plain matrices, no advanced structures). So you can read one element
    after the other linearily from the stream. Error checking is done
    per element.

    See <em>"Level 4 MAT-file format"</em> in the MATLAB documentation.

    ## Usage

    - Most errors are flagged by raising a VC::base::vc_runtime_error
      exception.

    - MatrixInfo::read() reads all information (header and name): you
      don't need to create MatrixHeader instances directly, but many
      methods for accessing attributes and reading data are
      implemented for MatrixHeader (and inherited by MatrixInfo). Use
      the various `expect_XXX` (e.g., MatrixInfo::expect_dimensions())
      calls for validation.

    - There are various functions in mat4 namespace for writing
      matrices, e.g., VC::math::mat4::write_matrix().

    - Directory reads information on matrix data from a MAT4 file and
      stores file offsets to data (for Directory::seekg() followed by,
      e.g., MatrixHeader::read_data() where Directory::Entry contains
      the header).  <em>Note: not for streaming.</em>

    - There are various stream output operators << based on, e.g.,
      MatrixInfo::to_s() methods.

    ## Notes

    - GNU Octave currently cannot read sparse matrices.
      We recommend, not to use write_sparse().

    - MATLAB seems to have trouble with data stored in big endian
      format (tested on little endian machine) wheras GNU Octave does
      fine.  (The documentation does not specify the endianess of the
      header, which encodes endianess!)

    - If you need to read sparse matrices into MATLAB, make sure
      indices are sorted correctly (e.g., by sort_sparse_ccs()).
      MATLAB's `load` expects coordinates in correct order, whereas
      write_sparse() does not sort data automatically.

    \sa VC::math::blas
 */

/** \example vclibs/math/example_mat4_read.cc
    \ingroup vc_mat4
    This is an example for how to use VC::math::mat4::MatrixInfo for reading.
    The program reads from stdin.
 */

/** \example vclibs/math/example_mat4_write.cc
    \ingroup vc_mat4
    This is an example for how to use VC::math::mat4 for writing.
    The program write some sample data to stdout.
 */

/** \example vclibs/math/example_mat4_directory.cc
    \ingroup vc_mat4
    This is an example for how to use VC::math::mat4::Directory.
 */

//=============================================================================

/// data format: type of matrix element \ingroup vc_mat4
enum DataFormat {
  Double=0,
  Single=1,
  Float=1,
  Int32=2,
  Int16=3,
  UInt16=4,
  UInt8=5,
  InvalidFormat=-1
};

/// matrix type \ingroup vc_mat4
enum MatrixType {
  Dense=0, //!< dense numeric matrix
  Text=1,  //!< text encoded as matrix
  Sparse=2 //!< sparse matrix encoded as m x 3 coordinate matrix
};


/// map DataFormat to C++ type \ingroup vc_mat4
template <DataFormat F> struct type_map {
# ifdef DOXYGEN_SKIP
  typedef NATIVE_TYPE value_type; //!< C++ type
# endif
};

# ifndef DOXYGEN_SKIP
template <> struct type_map<Double> { typedef double value_type; };
template <> struct type_map<Single> { typedef float value_type; };
template <> struct type_map<Int32>  { typedef int32_t value_type; };
template <> struct type_map<Int16>  { typedef int16_t value_type; };
template <> struct type_map<UInt16> { typedef uint16_t value_type; };
template <> struct type_map<UInt8>  { typedef uint8_t value_type; };
# endif

/// map type to DataFormat  \ingroup vc_mat4
template <typename T>
struct format_map {
  /// DataFormat
  static DataFormat id() { return InvalidFormat; }
};

# ifndef DOXYGEN_SKIP
template <> struct format_map<double>  { static DataFormat id() { return Double; } };
template <> struct format_map<float>   { static DataFormat id() { return Single; } };
template <> struct format_map<int32_t> { static DataFormat id() { return Int32; } };
template <> struct format_map<int16_t> { static DataFormat id() { return Int16; } };
template <> struct format_map<uint16_t>{ static DataFormat id() { return UInt16; } };
template <> struct format_map<uint8_t> { static DataFormat id() { return UInt8; } };
# endif

/// map DataFormat to VC::base::primtype:ID \ingroup vc_mat4
base::primtype::ID primtype_id(DataFormat _format);

/** \class MatrixHeader
    \brief MAT4 matrix header
    \ingroup vc_mat4

    **Note:** you may not want to use this class directly but
    prefer using MatrixInfo.

    \sa MatrixInfo
 */
class MatrixHeader {
public:
  /// invalid (empty) matrix
  MatrixHeader() { std::memset(data,0,sizeof(data)); }
  /// Copy constructor
  MatrixHeader(const MatrixHeader& _other) {
    assert(sizeof(data)==sizeof(MatrixHeader));
    std::memcpy(data,_other.data,sizeof(data));
  }

  /** Setup matrix header.
      Setup for write(), checks parameters.
      \param _byte_order byte order
      \param _format data format
      \param _type matrix type
      \param _mrows number of rows
      \param _ncols number of columns
      \param _imagf Has imaginary part?
      \param _namlen length of name + 1
      \throw VC::base::vc_runtime_error if !is_valid()
  */
  MatrixHeader(base::endian::ByteOrder _byte_order,
               DataFormat _format,
               MatrixType _type,
               int _mrows,int _ncols,bool _imagf,int _namlen)
    throw(VC::base::vc_runtime_error);


  /** @name Fields as defined in MATLAB documentation
      @{
  */

  int type_mopt() const { return data[0]; }
  int type_m() const { return data[0]/1000; }
  int type_o() const { return (data[0]/100)%10; }
  int type_p() const { return (data[0]/10)%10; }
  int type_t() const { return data[0]%10; }

  int mrows() const { return data[1]; }
  int ncols() const { return data[2]; }
  int imagf() const { return data[3]; }
  int namlen() const { return data[4]; }

  /// @}

  /** @name query fields (all assume is_valid())
      @{
  */

  /// get byte order (from type_m())
  base::endian::ByteOrder byte_order() const {
    return type_m()==0 ? base::endian::Little : base::endian::Big;
  }
  /// get data format (from type_p())
  DataFormat data_format() const { return DataFormat(type_p()); }

  /// get matrix type (from type_t())
  MatrixType matrix_type() const { return MatrixType(type_t()); }

  /// get textual description of _format ("double","float","int32_t",...)
  static const char* data_format_str(DataFormat _format);

  /// get textual description of _mtype ("dense","text","sparse")
  static const char* matrix_type_str(MatrixType _mtype);

  /// sparse matrix encoded as m x 3 or m x 4 matrix
  bool is_sparse() const { return type_t()==Sparse; }
  /// text matrix
  bool is_text() const { return type_t()==Text; }
  /// !is_complex()
  bool is_real() const {
    return (!is_sparse() && imagf()==0) || (is_sparse() && ncols()==3);
  }
  /// Has imaginary part? (imagf() or sparse matrix with 4 columns)
  bool is_complex() const { return !is_real(); }

  /// get string representation (for pretty-printed output)
  std::string to_s() const;

  /// @}

  /** @name check header
      @{
  */

  /// check header
  bool is_valid() const { return error()==0; }
  /// rerturns description error of !is_valid() or 0 otherwise
  const char* error() const;

  /// @}

  /** @name query size (all assume is_valid())
      @{
  */

  /// get size of a single data element (e.g., float for data_format()==Single)
  size_t data_element_size() const;

  /// number of real values (includes 4th=imaginary column for sparse)
  size_t n_real() const { return mrows()*ncols(); }

  /// number of imaginary values (always 0 for is_sparse()==true)
  size_t n_imag() const { return imagf()==0 ? 0 : mrows()*ncols(); }

  //// get number of data elements (double,float,...)
  size_t n_numbers() const {
    return n_real()+n_imag();
  }

  /// get data size in bytes (matrix data w/o namlen())
  size_t data_size() const { return n_numbers()*data_element_size(); }

  /// get total number of bytes (data_size()+namlen)
  size_t size() const { return data_size()+namlen(); }

  /// @}

  /** @name Input

      The following applies to all methods.

      - All methods throw VC::base::vc_runtime_error on failure.

      - The caller is responsible for allocating the memory buffers
        (e.g., `_real` and `_imag` arguments).

      - Some methods take an `_imag` argument for th imaginary part of
        a matrix. The following rule applies: If is_complex()==true
        then an `_imag` buffer **must** be provided, otherwise a
        vc_runtime_error is raised.  If the matrix is_real() then
        `_imag==0` is *ignored* and for `imag!=0` all (n_real())
        elements are set to *zero*.  This means if `_imag!=0` is
        provided, the caller \a must allocate a buffer regardless of
        is_complex()!

      - Sparse matrices are encoded as m x **N** or m x **N**
        matrices.  For real sparse matrices **N**=3, for complex
        sparse matrices **N**=4.  Consider using read_sparse().

      - is_complex() and is_real() respect sparseness: even though for
        sparse matrices we \a always have <tt>imagf()==0</tt>,
        is_complex() may be true.

      @{
  */

  /** Read from _in.
      At least sizeof(MatrixHeader) bytes are read from _in.
      On i/o error or if is_valid()==false after reading,
      an exception is thrown.
      \param _in input stream
      \return _in
      \b Note: you may not want to use this method directly but
      prefer MatrixInfo::read().
   */
  std::istream&
  read(std::istream& _in) throw(VC::base::vc_runtime_error);

  /** Read name.
      At least namlen() bytes are read from _in.
      On i/o error or if string was not null-terminated,
      an exception is thrown (and _name remains untouched).
      \param _in input stream
      \param[out] _name name of matrix
      \return _is
      \b Note: you may not want to use this method directly but
      prefer MatrixInfo::read().
  */
  std::istream& read_name(std::istream& _in,std::string& _name) const
    throw(VC::base::vc_runtime_error);

  /** Read data.
      Same as read_raw_data() but data elements will be converted to
      type T (without any checks).
      \sa read_matrix(), read_raw_data()
   */
  template <typename T>
  std::istream& read_data(std::istream& _in,T* _real,T* _imag=0) const
    throw(VC::base::vc_runtime_error) {
    check_real_imag(_real,_imag);
    read_some_data(_in,_real,n_real());
    if (imagf()==1)
      read_some_data(_in,_imag,n_imag());
    else if (_imag!=0) {
      for (int i=0;i<int(n_real());++i)
        _imag[i]=T(0);
    }
    return _in;
  }

  /** Read matrix.
      Same as read_data() with additional leading dimension parameter _ld.
      \param _in
      \param[out] _real real part as GE matrix (see, e.g., [\ref vc_blas])
      \param[out] _imag imaginary part as GE matrix
      \param _ld leading dimension of _real (and _imag), must be >= mrows()
      \return _is
   */
  template <typename T>
  std::istream& read_matrix(std::istream& _in,int _ld,T* _real,T* _imag=0) const
    throw(VC::base::vc_runtime_error);

  /** Read matrix and check dimensions.
      Same as read_matrix() but throw exception if matrix dimensions
      differ from _mrows x _ncols.
      \param _in
      \param _mrows expected # of rows
      \param _ncols expected # of columns
      \param[out] _real real part as GE matrix (see, e.g., [\ref vc_blas])
      \param[out] _imag imaginary part as GE matrix
      \param _ld leading dimension of _real (and _imag), must be >= mrows()
      \return _is
   */
  template <typename T>
  std::istream& read_matrix(std::istream& _in,int _mrows,int _ncols,int _ld,
                            T* _real,T* _imag=0) const
    throw(VC::base::vc_runtime_error);

  /// Read scalar _real as 1x1 matrix (read_matrix() with dimensions checked).
  template <typename T>
  std::istream& read_scalar(std::istream& _in,T& _real) const
    throw(VC::base::vc_runtime_error) {
    return read_matrix(_in,1,1,1,&_real);
  }
  /// Read scalar (_real,_imag) as 1x1 matrix (read_matrix() with dimensions checked).
  template <typename T>
  std::istream& read_scalar(std::istream& _in,T& _real,T& _imag) const
    throw(VC::base::vc_runtime_error) {
    return read_matrix(_in,1,1,1,&_real,&_imag);
  }

  /** Read text _data into string.
      Assumes lines are stores in rows.
      Assumes and checks is_text()==true. Checks range of matrix elements.
      \param _in input stream
      \param[out] _text text from matrix, rows are separated by _delimiter,
      trailing blanks are removed
      \param _delimiter separates rows in output
      \param _lines_in_columns indicates that text is stores column-wise
      \return _is
   */
  std::istream& read_text(std::istream& _in,std::string& _text,
                          const std::string& _delimiter="\n",
                          bool _lines_in_columns=false) const
    throw(VC::base::vc_runtime_error);

  /** Read sparse matrix.
      Assumes and check is_sparse()==true. Checks range of indices.
      Indices _i, _j are shifted to base zero.
      \param _in input stream
      \param[out] _i first coordinate
      \param[out] _j second coordinate
      \param[out] _m sparse matrix dimension (# of rows)
      \param[out] _n sparse matrix dimension (# of columns)
      \param[out] _real values
      \param[out] _imag imaginary values
      \return _in
   */
  template <typename T>
  std::istream& read_sparse(std::istream& _in,int* _i,int* _j,int& _m,int& _n,
                            T* _real,T* _imag=0) const
    throw(VC::base::vc_runtime_error);

  /// @}

  /** @name Output

      Generally, you will not use these methods directly but prefer
      functions
      - VC::math::mat4::write_matrix()
      - VC::math::mat4::write_text()
      - VC::math::mat4::write_sparse()

      @{
  */

  /** Setup write header \a and name to _out.
      \param _out output stream
      \param _name matrix name (must match namlen()!)
      \return _out
   */
  std::ostream& write(std::ostream& _out,const std::string& _name) const
    throw(VC::base::vc_runtime_error);


  /// Write buffers _real (and _imag) to _out.
  template <typename T>
  std::ostream& write_data(std::ostream& _out,
                           const T* _real,const T*_imag=0) const
    throw(VC::base::vc_runtime_error) {
    check_real_imag(_real,_imag);
    write_some_data(_out,_real,n_real());
    if (imagf()==1)
      write_some_data(_out,_imag,n_imag());
    return _out;
  }

  /// @}


  /** @name Low-level routines

      Generally, you will not need to call these methods directly.

      @{
  */

  /// Read _n data elements (helper called by read_data())
  template <typename T>
  std::istream& read_some_data(std::istream& _in,T* _buf,size_t _n) const
    throw(VC::base::vc_runtime_error) {
    return base::primtype::read(_in,_buf,_n,
                                primtype_id(data_format()),byte_order());
  }

  /// Write _n data elements (helper)
  std::ostream&
  write_some_data(std::ostream& _out,
                  DataFormat _buf_format,const void* _buf,size_t _n) const
    throw(VC::base::vc_runtime_error) {
    return base::primtype::write(_out,primtype_id(_buf_format),(const char*) _buf,_n,
                                 primtype_id(data_format()),byte_order());
  }

  /// Write _n data elements (helper called by write_data())
  template <typename T>
  std::ostream&
  write_some_data(std::ostream& _out,const T* _buf,size_t _n) const
    throw(VC::base::vc_runtime_error) {
    return base::primtype::write(_out,_buf,_n,
                                 primtype_id(data_format()),byte_order());
  }

  /// Read sparse coordinates as integers
  std::istream&
  read_sparse_coordinates(std::istream& _in,int* _i,int* _j,int& _m,int& _n) const
    throw(VC::base::vc_runtime_error);

  /// called by read_data(), write_data(), write_matrix(), throw
  void check_real_imag(const void* _real,const void* _imag) const
    throw(VC::base::vc_runtime_error);

  /// @}

private:

  uint32_t data[5];

  enum { BUFSIZE=64 /*!< buffer size for IO conversion */ };
};

//-----------------------------------------------------------------------------

/** \class MatrixInfo
    \brief Information about MAT4 matrix.
    \ingroup vc_mat4
    Provides MatrixHeader and name() (from MatrixHeader::read()).
    This is all you need to read matrix data (using
    MatrixHeader::read_data(),MatrixHeader::read_matrix(),...).
    For output use functions VC::math::mat4::write_matrix(),
    VC::math::mat4::write_text(),VC::math::mat4::write_sparse().
    \sa MatrixInfo
 */
class MatrixInfo : public MatrixHeader {
public:
  /// empty
  MatrixInfo() {}
  /// Copy constructor
  MatrixInfo(const MatrixInfo& _other)
  : MatrixHeader(_other), m_name(_other.m_name) {}

  /// calls read()
  MatrixInfo(std::istream& _in) throw(VC::base::vc_runtime_error) {
    read(_in);
  }
  /// setup from arguments (see MatrixHeader::MatrixHeader())
  MatrixInfo(const std::string& _name,
             int _mrows,int _ncols,
             MatrixType _type,DataFormat _format,
             bool _imgaf=false,
             VC::base::endian::ByteOrder _byte_order=VC::base::endian::Native)
  throw(VC::base::vc_runtime_error)
    : MatrixHeader(_byte_order,_format,_type,
                   _mrows,_ncols,_imgaf,int32_t(_name.length())+1),
      m_name(_name) {}

  /// get name
  const std::string& name() const { return m_name; }

  /// get string representation (for pretty-printed output)
  std::string to_s() const {
    return m_name+' '+MatrixHeader::to_s();
  }

  /** @name Input

      \b Note: use methods inherited from MatrixHeader for reading
      data, e.g., MatrixHeader::read_data().

      @{
  */

  /** Read from _in.
      Calls MatrixHeader::read(): reads header and name.
      \param _in input stream
      \return _in
      \sa MatrixHeader::read(), MatrixHeader::read_name().
   */
  std::istream&
  read(std::istream& _in) throw(VC::base::vc_runtime_error) {
    MatrixHeader::read(_in);
    return read_name(_in,m_name);
  }

  /** read() matrix headers until matrix _name is found.
      \param _in input stream (read() headers,
      data is skipped using istream::seekg()
      \param _name of data
      \return _in
      \throw vc_runtime_error on error or if _in.eof()
   */
  std::istream&
  skip_until(std::istream& _in,
            const std::string& _name) throw(VC::base::vc_runtime_error);

  /// @}

  /** @name Validation

      Throw assertion if expectations are not met. Use for validation
      after read().

      @{
   */

  /// expect _mtype
  void expect_matrix_type(MatrixType _type) const
    throw(VC::base::vc_runtime_error);
  /// expect Dense matrix
  void expect_dense() const throw(VC::base::vc_runtime_error) {
    expect_matrix_type(Dense);
  }
  /// expect Text matrix
  void expect_text() const throw(VC::base::vc_runtime_error) {
    expect_matrix_type(Text);
  }
  /// expect Sparse matrix
  void expect_sparse() const throw(VC::base::vc_runtime_error) {
    expect_matrix_type(Sparse);
  }
  /// expect real matrix
  void expect_real(bool _real=true) const
    throw(VC::base::vc_runtime_error);
  /// expect complex matrix (with mandatory imaginary part)
  void expect_complex(bool _complex=true) const
    throw(VC::base::vc_runtime_error);

  /** Expect matrix dimensions.
      Not for sparse matrices, their dimension is available only after reading
      the data by read_sparse().
   */
  void expect_dimensions(int _m,int _n) const
    throw(VC::base::vc_runtime_error);

  /// expect vector (mrows()==1 || nrcols()==1)
  void expect_vector(int _n=-1) const throw(VC::base::vc_runtime_error);

  /// expect scalar
  void expect_scalar() const throw(VC::base::vc_runtime_error);

  /// expect _m mrows()
  void expect_mrows(int _m) const throw(VC::base::vc_runtime_error);

  /// expect _n ncols()
  void expect_ncols(int _n) const throw(VC::base::vc_runtime_error);

  /// @}

private:
  std::string      m_name; //!< name
};

//-----------------------------------------------------------------------------

/** \class Directory
    \brief Directory of a MAT4-file.
    \ingroup vc_mat4

    Use this class to read all MatrixInfo from a MAT4-file (std::istream) and
    storing the file offsets of the data parts (Directory::Entry).

    The data is \a not read, instead std::istream::seekg() is used to
    get to the next MatrixInfo. The idea is to first check for availability
    and correct type of matrices and to then read the data parts.

    ##Notes

    - If two entries share the same name then only the last entry is
      stored in map(). The default setting is to raise an exception,
      see Directory::read(). (Note that index() then contains multiple
      instances of the same iterator on map().)
    - Reading the directory and then the data would not work for in a
      streaming manner, e.g., for `cin`. The input stream must
      support seeking in both directions.

    \sa MatrixInfo
 */
class Directory {
public:

  /// Directory entry (MatrixInfo and data_pos) \ingroup vc_mat4
  struct Entry : public MatrixInfo {
    /// default constructor: invalid
    Entry() : data_pos(0) {}
    /// copy constructor
    Entry(const Entry& _other)
      : MatrixInfo(_other), data_pos(_other.data_pos) {}
    /// construct from _in: call MatrixInfo::read(), set data_pos
    Entry(std::istream& _in) throw(VC::base::vc_runtime_error);

    /// advance get pointer of _in to data_pos
    std::istream& seekg(std::istream& _in) const
      throw(VC::base::vc_runtime_error);

    /// advance get pointer of _in to header_pos
    std::istream& seekg_header(std::istream& _in) const
      throw(VC::base::vc_runtime_error);

    /** Get string representation (for pretty-printed output).
        Appends data_pos as hex number.
     */
    std::string to_s() const;

    /// get stream position of MatrixHeader from data_pos
    std::streampos header_pos() const {
      return data_pos-std::streamoff(namlen()+sizeof(MatrixHeader));
    }

    /// stream position of data (after MatrixInfo stored in this instance)
    std::streampos data_pos;
  };

  /// map() type
  typedef std::map<std::string,Entry> map_type;
  /// iterator on map()
  typedef map_type::const_iterator const_iterator;

  /// index() type (stores const_iterator in order of appearance)
  typedef std::vector<const_iterator> index_type;
  /// interator on index()
  typedef index_type::const_iterator index_iterator;

  /// default constructor (empty directory)
  Directory() : m_in(0), m_base(0), m_end(0) {}

  /** Construct from _in.
      Calls read(). _in is used only for read(), it is not closed on
      destruction.
  */
  Directory(std::istream& _in,bool _check_unique=true)
    throw(VC::base::vc_runtime_error);

  /// get attached input stream
  std::istream& input() const throw(VC::base::vc_runtime_error) {
    if (m_in==0)
      throw VC_RUNTIME_ERROR("no input stream");
    return *m_in;
  }
  /// attach other std::istream (e.g., after a re-open); use with care
  std::istream* replace_input(std::istream& _new_in) {
    std::istream* old_in=m_in;
    m_in=&_new_in;
    return old_in;
  }

  /// get map_type instance which stores the mapping
  const map_type& map() const { return m_map; };

  /** Get file order index.

      index() returns an array (index_type) with const_iterator objects
      pointing to Enry instances in map(). In contrast to map(), the
      order of elements in index() is based on ascending Entry::data_pos.

      \b Note: see also Note on _check_unique argument in read()!

      \return reference to index_type object
  */
  const index_type& index() const { return m_index; };
  /// same as map()
  const map_type& operator()() const { return m_map; }

  /// get Entry _name (or throw exception, map().find() does not throw)
  const Entry& operator[](const std::string& _name) const
    throw(VC::base::vc_runtime_error);

  /// advance get pointer of in() to \a data of _ii (throw if _ii==map().end())
  std::istream& seekg(const_iterator _ii) const
    throw(VC::base::vc_runtime_error) {
    if (_ii==m_map.end())
      throw VC_RUNTIME_ERROR("invalid entry");
    return _ii->second.seekg(*m_in);
  }
  /// advance get pointer of in() to \a data _name
  std::istream& seekg(const std::string& _name) const
    throw(VC::base::vc_runtime_error) {
    return seekg(m_map.find(_name));
  }

  /// advance get pointer of in() to \a header of _ii (throw if _ii==map().end()
  std::istream& seekg_header(const_iterator _ii) const
    throw(VC::base::vc_runtime_error) {
    if (_ii==m_map.end())
      throw VC_RUNTIME_ERROR("invalid entry");
    return _ii->second.seekg_header(*m_in);

  }
  /// advance get pointer of in() to \a header _name
  std::istream& seekg_header(const std::string& _name) const
    throw(VC::base::vc_runtime_error) {
    return seekg_header(m_map.find(_name));
  }

  /// get stream position at call to read() (or construction)
  std::streampos base_pos() const { return m_base; }
  /// get stream position where read() stopped
  std::streampos end_pos() const { return m_end; }
  /// get size of data end_pos()-base_pos()
  std::streamoff data_size() const { return m_end-m_base; }

  /** Read data from _in.
      Attaches _in so it remains accessible via in(). Reads all MatrixInfo
      from _in and sets (istream::seekg()) get pointer back to base_pos().
      \param _in input stream
      \param _check_unique throw VC::base::vc_runtime_error if names are
      not unique, i.e., if a name appears twice.
      \return _in
      **Note:** If _check_unique==false and a name appears multiple times,
      then
      - only the *last* entry is stored in map()!
      - index() contains multiple instances of the same iterators,
        i.e., `index().size()>map().size()`, which all point to this
        *last* entry!
   */
  std::istream& read(std::istream& _in,bool _check_unique=true)
    throw(VC::base::vc_runtime_error);

private:
  Directory(const Directory&) {} //!< Directory is not copy-able!

  map_type       m_map;   //!< map name to MatrixInfo
  index_type     m_index; //!< stores order of MatrixInfo objects from read()

  std::istream*  m_in;    //!< current input stream
  std::streampos m_base;  //!< stream position where read() started
  std::streampos m_end;   //!< stream position where read() stopped
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/// print MatrixInfo::to_s() \ingroup vc_mat4
std::ostream& operator<<(std::ostream& _out,const MatrixInfo& _m);

/// print MatrixHeader::to_s() \ingroup vc_mat4
std::ostream& operator<<(std::ostream& _out,const MatrixHeader& _m);

/// print Directory::Entry::to_s() \ingroup vc_mat4
std::ostream& operator<<(std::ostream& _out,const Directory::Entry& _m);

/** print Directory _dir.
    \ingroup vc_mat4
    Print all entries as Directory::Entry::to_s(), each on a new line
    \param _out output
    \param _dir directory
    \return _out
*/
std::ostream& operator<<(std::ostream& _out,const Directory& _dir);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Write dense numeric matrix.
    \ingroup vc_mat4
    \param _out output stream
    \param _name matrix name
    \param _mrows number of rows
    \param _ncols number of columns
    \param _ld leading dimension of _real and _imag
    \param _src_format data format for _real, _imag
    \param _real real part of matrix
    \param _imag imaginary part of matrix, real matrix if _imag==0
    \param _dst_format data format written to _out
    \param _byte_order byte order written to _out
 */
std::ostream&
write_matrix(std::ostream& _out,const std::string& _name,
             int _mrows,int _ncols,int _ld,DataFormat _src_format,
             const void* _real,const void* _imag,DataFormat _dst_format,
             VC::base::endian::ByteOrder _byte_order=VC::base::endian::Native)
  throw(VC::base::vc_runtime_error);

/// Same as above for real matrix (_imag=0). \ingroup vc_mat4
inline std::ostream&
write_matrix(std::ostream& _out,const std::string& _name,int _mrows,int _ncols,int _ld,
             DataFormat _src_format,const void* _real,DataFormat _dst_format,
             VC::base::endian::ByteOrder _byte_order=VC::base::endian::Native)
  throw(VC::base::vc_runtime_error) {
  return write_matrix(_out,_name,_mrows,_ncols,_ld,
                      _src_format,_real,(const void*) 0,_dst_format,_byte_order);
}

/** Write dense numeric matrix.
    \ingroup vc_mat4
    \param _out output stream
    \param _name matrix name
    \param _mrows number of rows
    \param _ncols number of columns
    \param _ld leading dimension of _real and _imag buffers
    \param _real real part of matrix
    \param _imag imaginary part of matrix, real matrix if _imag==0
    \param _byte_order byte order written to _out
*/
template <typename T>
std::ostream&
write_matrix(std::ostream& _out,const std::string& _name,int _mrows,int _ncols,int _ld,
             const T* _real,const T* _imag=0,
             VC::base::endian::ByteOrder _byte_order=VC::base::endian::Native)
  throw(VC::base::vc_runtime_error) {
  return write_matrix(_out,_name,_mrows,_ncols,_ld,
                      format_map<T>::id(),_real,_imag,format_map<T>::id(),_byte_order);
}

/// Write scalar _real as 1x1 matrix (write_matrix()). \ingroup vc_mat4
template <typename T>
std::ostream&
write_scalar(std::ostream& _out,const std::string& _name,const T& _real) {
  T re=_real;
  return write_matrix(_out,_name,1,1,1,&re);
}

/// Write scalar (_real,_imag) as 1x1 matrix (write_matrix()). \ingroup vc_mat4
template <typename T>
std::ostream&
write_scalar(std::ostream& _out,
             const std::string& _name,const T& _real,const T& _imag) {
  T re=_real, im=_imag;
  return write_matrix(_out,_name,1,1,1,&re,&im);
}

/** Write _text.
    \ingroup vc_mat4
    \param _out output stream
    \param _name matrix name
    \param _text text data written as a single matrix row
    \return _out
   */
std::ostream&
write_text(std::ostream& _out,const std::string& _name,
           const std::string& _text)
  throw(VC::base::vc_runtime_error);



# ifndef DOXYGEN_SKIP
 /// comparison operator for CCS sort (see sort_sparse_ccs()) \ingroup vc_mat4
template <typename IDX>
struct CompareIdx_CCS {
  CompareIdx_CCS(const IDX* _i,const IDX* _j) : pi(_i), pj(_j) {}
  const IDX* pi;
  const IDX* pj;
  bool operator()(int _ka,int _kb) const {
    if      (pj[_ka]< pj[_kb]) return true;
    else if (pj[_ka]!=pj[_kb]) return false;

    if      (pi[_ka]< pi[_kb]) return true;
    else                       return false;
  }
};
# endif

/** Generate permutation for sorting sparse indices.
    \ingroup vc_mat4
    Called by sort_sparse_ccs().
    \param _n number of indeices (_i,_j)
    \param _i array of row indices
    \param _j array of column indices
    \param[out] _permutation permutation of _i, _j
    and as well real and imaginary parts
    \sa sort_sparse_ccs()
 */
template <typename IDX>
void
compute_ccs_permutation(int _n,const IDX* _i,const IDX* _j,int* _permutation) {
  for (int i=0;i<_n;++i)
    _permutation[i]=i;

  CompareIdx_CCS<IDX> cmp(_i,_j);
  std::sort(_permutation+0,_permutation+_n,cmp);
}

/** Apply permutation _p to (_i,_j,_real,_imag).
    \ingroup vc_mat4
    \param _n number of elements in _p and (_i,_j,_real,_imag)
    \param _p permitation
    \param _i row indices
    \param _j column indices
    \param _real real part
    \param _imag imaginary part, may be **0**
    \param _buf temporary buffer must equal **0** or provide size
    `_n*max(sizeof(IDX),sizeof(T))`
    \sa compute_ccs_permutation(), sort_sparse_ccs()
 */
template <typename IDX,typename T>
void
apply_permutation(int _n,const int* _p,IDX* _i,IDX* _j,T* _real,T* _imag=0,void* _buf=0) {
  const int*   p=_p;

  void*        malloced_buf=0;
  size_t       sz=_n*std::max(sizeof(IDX),sizeof(T));

# if VC_ALLOCA_LIMIT>0
  if (_buf==0 && sz<=VC_ALLOCA_LIMIT)
    _buf=alloca(sz);
# endif

  if (_buf==0)
    _buf=malloced_buf=std::malloc(sz);

  // Implementation note:
  //  using the extra memory buffer is significantly faster than swaps on p and reverse(p).

  memcpy(_buf,_i,_n*sizeof(IDX));
  for (int i=0;i<_n;++i)
    _i[i]=((IDX*) _buf)[p[i]];

  memcpy(_buf,_j,_n*sizeof(IDX));
  for (int i=0;i<_n;++i)
    _j[i]=((IDX*) _buf)[p[i]];

  memcpy(_buf,_real,_n*sizeof(T));
  for (int i=0;i<_n;++i)
    _real[i]=((T*) _buf)[p[i]];

  if (_imag!=0) {
    memcpy(_buf,_imag,_n*sizeof(T));
    for (int i=0;i<_n;++i)
      _imag[i]=((T*) _buf)[p[i]];
  }

  if (malloced_buf!=0)
    std::free(malloced_buf);
}

/** Sort sparse entries (_i,_j,_real,_imag) CCS-wise.
    \ingroup vc_mat4
    Enables construction of compressed column storage (CCS) sparse
    matrix. MATLAB's \c load expects sparse data to be sorted.
    \param _n number of entries (generally nnz)
    \param _i row index
    \param _j column index
    \param _real real data
    \param _imag imaginary data (may be **0**)
    \param _buf temporary buffer must equal **0** or provide size
    `_n*(sizeof(int)+max(sizeof(IDX),sizeof(T)))`
    \sa compute_ccs_permutation(), write_sparse()
 */
template <typename IDX,typename T>
void sort_sparse_ccs(int _n,IDX* _i,IDX* _j,T* _real,T* _imag=0,void* _buf=0) {
  void*        malloced_buf=0;
  size_t       sz=_n*(sizeof(int)+std::max(sizeof(IDX),sizeof(T)));

# if VC_ALLOCA_LIMIT>0
  if (_buf==0 && sz<=VC_ALLOCA_LIMIT)
    _buf=alloca(sz);
# endif

  if (_buf==0)
    _buf=malloced_buf=malloc(sz);

  int* p=(int*) _buf;

  // compute permutation p
  compute_ccs_permutation(_n,_i,_j,p);
  // apply permutation
  apply_permutation(_n,p,_i,_j,_real,_imag,(void*) (((int*) _buf)+_n));

  if (malloced_buf!=0)
    free(malloced_buf);
}

/** Write sparse matrix.
    \ingroup vc_mat4
    The matrix is \a always written in Double format.

    \b Note: MATLAB's \c load requires <tt>(i,j,re,im)</tt> sorted by index
    (CCS). If indices are not sorted, you must call sort_sparse_ccs() priort to
    write_sparse()!

    \param _out output stream
    \param _name matrix name
    \param _m # of rows in sparse matrix
    \param _n # of columns in sparse matrix
    \param _nnz number of enries in (_i,_j,_real,_imag)
    \param _i first coordinate (0 based)
    \param _j second coordinate (0 based)
    \param _real real part of matrix
    \param _imag imaginary part of matrix, real matrix if _imag==0
    \param _byte_order byte order written to _out
    \sa sort_sparse_ccs()
 */
template <typename IDX,typename T>
std::ostream&
write_sparse(std::ostream& _out,const std::string& _name,int _nnz,int _m,int _n,
             const IDX* _i,const IDX* _j,const T* _real,const T* _imag=0,
             VC::base::endian::ByteOrder _byte_order=VC::base::endian::Native)
  throw(VC::base::vc_runtime_error) {
  return write_sparse(_out,_name,_nnz,_m,_n,
                      format_map<IDX>::id(),_i,_j,
                      format_map<T>::id(),_real,_imag,_byte_order);
}

/** Write sparse matrix.
    \ingroup vc_mat4
    The matrix is \a always written in Double format.

    \b Note: MATLAB's \c load requires <tt>(i,j,re,im)</tt> sorted by index
    (CCS). If indices are not sorted, you must call sort_sparse_ccs() priort to
    write_sparse()!

    \param _out output stream
    \param _name matrix name
    \param _m # of rows in sparse matrix
    \param _n # of columns in sparse matrix
    \param _nnz number of enries in (_i,_j,_real,_imag)
    \param _idx_format format of index buffers _i,_j
    \param _i first coordinate (0 based)
    \param _j second coordinate (0 based)
    \param _data_format format of  buffers _real,_imag
    \param _real real part of matrix
    \param _imag imaginary part of matrix, real matrix if _imag==0
    \param _byte_order byte order written to _out
    \sa sort_sparse_ccs()
 */
std::ostream&
write_sparse(std::ostream& _out,const std::string& _name,int _m,int _n,int _nnz,
             DataFormat _idx_format,const void* _i,const void* _j,
             DataFormat _data_format,const void* _real,const void* _imag=0,
             VC::base::endian::ByteOrder _byte_order=VC::base::endian::Native)
  throw(VC::base::vc_runtime_error);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template <typename T>
std::istream&
MatrixHeader::read_matrix(std::istream& _in,int _ld,T* _real,T* _imag) const
  throw(VC::base::vc_runtime_error) {

  if (_ld<mrows())
    throw VC_RUNTIME_ERROR("invalid leading dimension");

  if (_ld==mrows())
    return read_data(_in,_real,_imag);

  check_real_imag(_real,_imag);

  for (int j=0;j<ncols();++j)
    read_some_data(_in,_real+j*_ld,mrows());

  if (imagf()==1) {
    for (int j=0;j<ncols();++j)
      read_some_data(_in,_imag+j*_ld,mrows());
  }
  else if (_imag!=0) {
    for (int j=0;j<ncols();++j)
      for (int i=0;i<mrows();++i)
        _imag[j*_ld+i]=T(0);
  }

  return _in;
}
//-----------------------------------------------------------------------------

template <typename T>
std::istream&
MatrixHeader::read_matrix(std::istream& _in,
                          int _mrows,int _ncols,int _ld,T* _real,T* _imag) const
  throw(VC::base::vc_runtime_error) {

  if (mrows()!=_mrows || ncols()!=_ncols) {
    std::string msg=VC::base::formatf("expected %d x %d matrix but got %d x %d",
                                      _mrows,_ncols,mrows(),ncols());
    throw VC_RUNTIME_ERROR(msg);
  }
  return read_matrix(_in,_ld,_real,_imag);
}

//-----------------------------------------------------------------------------

template <typename T>
std::istream&
MatrixHeader::read_sparse(std::istream& _in,int* _i,int* _j,int& _m,int& _n,
                          T* _real,T* _imag) const
  throw(VC::base::vc_runtime_error) {

  check_real_imag(_real,_imag);
  read_sparse_coordinates(_in,_i,_j,_m,_n);
  read_some_data(_in,_real,mrows());
  if (ncols()==4) {
    assert(_imag!=0 && "checked above");
    read_some_data(_in,_imag,mrows());
  }
  else if (_imag!=0) {
    for (int i=0;i<mrows();++i)
      _imag[i]=T(0);
  }
  return _in;
}

//-----------------------------------------------------------------------------

# if defined(_MSC_VER)
#  pragma warning (pop)
# endif

//=============================================================================
} // namespace mat4io
} // namespace base
} // namespace VC

# endif // VC_BASE_MAT4_HH
