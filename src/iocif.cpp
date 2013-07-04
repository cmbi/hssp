#include "mas.h"

#include <cassert>

#include <iostream>
#include <string>
#include <vector>

#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "utils.h"
#include "iocif.h"

using namespace std;
namespace io = boost::iostreams;

//	Our CIF implementation consists of flyweight classes.

namespace
{

// skip routines to quickly position character pointer p at a next interesting location
const char* skip_line(const char* p, const char* end);
const char* skip_white(const char* p, const char* end);	// skip over white-space and comments
const char* skip_value(const char* p, const char* end);

enum CIF_file_type { cif_char, cif_number, cif_text };

struct CIF_field;
struct CIF_record;

//struct CIF_value
//{
//	CIF_record*		m_record;
//	CIF_field*		m_field;
//
//	string			str() const;
//	double			number() const;
//	
//};
//
struct CIF_field
{
	string			name() const			{ return string(m_name, m_name_end); }
	string			value() const			{ return string(m_data, m_data_end); }
	
	bool operator==(const CIF_field& rhs) const	{ return m_name == rhs.m_name and m_data == rhs.m_data; }

	const char*		m_name;
	const char*		m_name_end;
	const char*		m_data;
	const char*		m_data_end;
	CIF_file_type	m_type;
};

struct CIF_row
{
	CIF_field operator[](const char* inName) const;
	bool operator==(const CIF_row& rhs) const		{ return m_fields == rhs.m_fields; }
	
	vector<CIF_field>	m_fields;
};

CIF_field CIF_row::operator[](const char* inName) const
{
	foreach (const CIF_field& f, m_fields)
	{
		if (strncmp(inName, f.m_name, f.m_name_end - f.m_name) == 0)
			return f;
	}
	
	throw mas_exception(boost::format("Field %s not found") % inName);
	return CIF_field();
}

struct CIF_record
{
	string			name() const			{ return m_name; }
	
	struct const_iterator : public iterator<forward_iterator_tag,const CIF_row>
	{
		typedef iterator<forward_iterator_tag, const CIF_row>	base_type;
		typedef base_type::reference							reference;
		typedef base_type::pointer								pointer;
		
						const_iterator();
						const_iterator(const CIF_record& rec, const CIF_row& row) : m_rec(rec), m_row(row) {}
						const_iterator(const const_iterator& iter) : m_rec(iter.m_rec), m_row(iter.m_row) {}
		const_iterator&	operator=(const const_iterator& iter)			{ m_row = iter.m_row; return *this; }

		reference		operator*() const								{ return m_row; }
		pointer			operator->() const								{ return &m_row; }

		const_iterator&	operator++()									{ m_rec.advance(m_row); return *this; }
		const_iterator	operator++(int)									{ const_iterator iter(*this); operator++(); return iter; }

		bool			operator==(const const_iterator& iter) const	{ return m_row == iter.m_row; }
		bool			operator!=(const const_iterator& iter) const	{ return not operator==(iter); }

	  private:
		const CIF_record&	m_rec;
		CIF_row				m_row;
	};
	
	CIF_row			front() const;
	CIF_row			back() const;
	
	const_iterator	begin() const;
	const_iterator	end() const;

	typedef const_iterator iterator;

	void			advance(CIF_row& row) const;	// update pointers to next data row, if any

	bool operator<(const CIF_record& rhs) const							{ return m_name < rhs.m_name; }
	
	const char*		m_start;
	const char*		m_end;
	bool			m_loop;
	uint32			m_field_count;
	string			m_name;
};

CIF_row CIF_record::front() const
{
	CIF_row result;
	
	const char* p = m_start;
	
	if (m_loop)
	{
		for (uint32 i = 0; i < m_field_count; ++i)
		{
			assert(*p == '_');
			assert(*(p + m_name.length()) == '.');
			
			CIF_field field = {};

			field.m_name = p = p + m_name.length() + 1;
			while (p != m_end and not isspace(*p))
				++p;
			
			field.m_name_end = p;
			
			p = skip_white(p, m_end);
			
			result.m_fields.push_back(field);
		}

		foreach (CIF_field& fld, result.m_fields)
		{
			fld.m_data = skip_white(p, m_end);
			fld.m_data_end = skip_value(fld.m_data, m_end);
			p = skip_white(fld.m_data_end, m_end);
		}
	}
	else
	{
		for (uint32 i = 0; i < m_field_count; ++i)
		{
			assert(*p == '_');
			assert(*(p + m_name.length()) == '.');
			
			CIF_field field;
			field.m_name = p = p + m_name.length() + 1;
			while (p != m_end and not isspace(*p))
				++p;
			
			field.m_name_end = p;
			
			p = skip_white(p, m_end);
			field.m_data = p;

			p = skip_value(p, m_end);
			field.m_data_end = p;
			
			result.m_fields.push_back(field);
		}	
	}
	
	return result;
}

CIF_record::iterator CIF_record::begin() const
{
	return const_iterator(*this, front());
}

CIF_record::iterator CIF_record::end() const
{
	return const_iterator(*this, CIF_row());
}

void CIF_record::advance(CIF_row& row) const
{
	if (m_loop and not row.m_fields.empty())
	{
		const char* p = skip_white(row.m_fields.back().m_data_end, m_end);

		if (p >= m_end)
			row.m_fields.clear();
		else
		{
			foreach (CIF_field& fld, row.m_fields)
			{
				fld.m_data = skip_white(p, m_end);
				fld.m_data_end = skip_value(fld.m_data, m_end);
				p = skip_white(fld.m_data_end, m_end);
			}
		}
	}
	else
		row.m_fields.clear();
}

struct CIF_file
{
	CIF_file(const char* inData, size_t inSize);

	CIF_record operator[](const char* inName) const;

	// skip to first character after the next NL character
	const char* skip_line(const char* p)		{ return ::skip_line(p, m_end); }
	
	// skip over values for a record
	const char* skip_value(const char* p)		{ return ::skip_value(p, m_end); }

	vector<CIF_record>	m_records;
	const char*			m_data;
	const char*			m_end;
};

CIF_file::CIF_file(const char* inData, size_t inSize)
	: m_data(inData), m_end(inData + inSize)
{
	// CIF files are simple to parse
	
	const char* p = m_data;
	
	if (strncmp(p, "data_", 5) != 0)
		throw mas_exception("Is this an mmCIF file?");
	
	p = skip_line(p);
	
	CIF_record rec = { p };
	uint32 valueCount = 0;
	bool loop = false;
	
	while (p < m_end)
	{
		if (isspace(*p))	// skip over white space
		{
			++p;
			continue;
		}
		
		if (*p == '#')	 // line starting with hash, this is a comment, skip
		{
			p = skip_line(p);
			continue;	
		}
		
		if (strncmp(p, "loop_", 5) == 0)
		{
			if (not m_records.empty() and m_records.back().m_end == nullptr)
				m_records.back().m_end = p;

			loop = true;
			rec.m_loop = false;
			p = skip_line(p + 5);

			continue;
		}
		
		const char* s = p;
		
		if (*p == '_')	// a label
		{
			// scan for first dot
			bool newName = loop;
			const char* n = rec.m_start;
			
			for (;;)
			{
				if (not newName and *p != *n)
					newName = true;
				
				++p;
				++n;
				
				if (p == m_end or *p == '.' or isspace(*p))
					break;
			}
			
			if (*p == '.')	// OK, found a record
			{
				if (newName)
				{
					// store start as end for the previous record, if any
					if (not m_records.empty() and m_records.back().m_end == nullptr)
						m_records.back().m_end = s;

					rec.m_start = s;
					rec.m_end = nullptr;
					rec.m_loop = loop;
					rec.m_field_count = 1;
					rec.m_name = string(s, p);

					m_records.push_back(rec);
				}
				else
					m_records.back().m_field_count += 1;
				
				// skip over field name
				while (p != m_end and not isspace(*p))
					++p;
			}
			else
			{
				// store start as end for the previous record, if any
				if (not m_records.empty() and m_records.back().m_end == nullptr)
					m_records.back().m_end = s;

				// a record without a field (is that possible in mmCIF?)
				cerr << "record without field: " << string(s, p) << endl;
				
				rec.m_start = s;
				rec.m_end = nullptr;
				rec.m_loop = loop;
				rec.m_field_count = 0;
				rec.m_name = string(s, p);

				m_records.push_back(rec);
			}

			if (not rec.m_loop)
				p = skip_value(p);
			
			loop = false;
			continue;
		}
		
		if (rec.m_loop == false)
		{
			// guess we should never reach this point
			throw mas_exception("invalid CIF file? (unexpected data, not in loop)");
		}
		
		p = skip_value(p);
		
		// check for a new data_ block
		if (p != m_end and strncmp(p, "data_", 5) == 0)
			throw mas_exception("Multiple data blocks in CIF file");
	}

	if (not m_records.empty() and m_records.back().m_end == nullptr)
		m_records.back().m_end = p;
	
	sort(m_records.begin(), m_records.end());
}

CIF_record CIF_file::operator[](const char* inName) const
{
	CIF_record test;
	test.m_name = inName;
	
	vector<CIF_record>::const_iterator i = lower_bound(m_records.begin(), m_records.end(), test);
	if (i == m_records.end())
		throw mas_exception(boost::format("Field %s not found") % inName);
	
	return *i;
}

// skip to first character after the next NL character
const char* skip_line(const char* p, const char* end)
{
	while (p != end)
	{
		if (*p++ == '\n')
			break;
	}

	return p;
}

// skip over white space and comments
const char* skip_white(const char* p, const char* end)
{
	while (p != end)
	{
		if (isspace(*p))
		{
			++p;
			continue;
		}
		
		if (*p == '#')
		{
			do ++p; while (p < end and *p != '\n');
			continue;
		}
		
		break;
	}

	return p;
}

// skip over values for a record
const char* skip_value(const char* p, const char* end)
{
	for (;;)
	{
		if (*p == '\n' and *(p + 1) == ';')
		{
			do p = skip_line(p + 1, end); while (p < end and *p != ';');
			++p;
			break;
		}
		
		if (isspace(*p))
		{
			++p;
			continue;
		}

		if (*p == '\'')
		{
			do ++p; while (p != end and *p != '\'');
			++p;
			break;
		}
		
		if (*p == '\"')
		{
			do ++p; while (p != end and *p != '\"');
			++p;
			break;
		}
		
		if (*p == '#')
		{
			p = skip_line(p + 1, end);
			continue;
		}
		
		while (p != end and not isspace(*p))
			++p;
		
		break;
	}
	
	return p;
}

	
}

void ReadCIF(std::istream& in, MProtein& out)
{
	vector<char> buffer;
	buffer.reserve(10 * 1024 * 1024);	// reserve 10 MB, should be sufficient for most

	io::copy(in, io::back_inserter(buffer));
	buffer.push_back(0);				// end with a null character, makes coding a bit easier
	
	CIF_file cif(&buffer[0], buffer.size() - 1);
	
//	foreach (const CIF_record& r, cif.m_records)
//		cout << r.name() << '\t' << r.m_field_count << endl;
	
	cout << "id: " << cif["_entry"].front()["id"].value() << endl;
	
	foreach (const CIF_row& row, cif["_atom_type"])
	{
		cout << row["symbol"].value() << endl;
	}

	foreach (const CIF_row& row, cif["_atom_site"])
	{
		cout << "ATOM  " << row["Cartn_x"].value() << ' ' << row["Cartn_y"].value() << ' ' << row["Cartn_z"].value() << endl;
	}
}

