/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class Checkpoint
{
private:
	const std::string _filename;
	FILE * const _file;
	uint32_t _crc32 = 0;

	// Rosetta Code, CRC-32, C
	static uint32_t rc_crc32(const uint32_t crc32, const char * const buf, const size_t len)
	{
		static uint32_t table[256];
		static bool have_table = false;

		// This check is not thread safe; there is no mutex
		if (!have_table)
		{
			// Calculate CRC table
			for (size_t i = 0; i < 256; ++i)
			{
				uint32_t rem = static_cast<uint32_t>(i);  // remainder from polynomial division
				for (size_t j = 0; j < 8; ++j)
				{
					if (rem & 1)
					{
						rem >>= 1;
						rem ^= 0xedb88320;
					}
					else rem >>= 1;
				}
				table[i] = rem;
			}
			have_table = true;
		}

		uint32_t crc = ~crc32;
		for (size_t i = 0; i < len; ++i)
		{
			const uint8_t octet = static_cast<uint8_t>(buf[i]);  // Cast to unsigned octet
			crc = (crc >> 8) ^ table[(crc & 0xff) ^ octet];
		}
		return ~crc;
	}

	void error(const std::string & str) const { std::cerr << "Error: file '" << _filename << "', " << str << "." << std::endl; }

public:
	Checkpoint(const std::string & filename, const char * const mode) : _filename(filename), _file(std::fopen(filename.c_str(), mode)), _crc32(0) {}
	virtual ~Checkpoint() { if (_file != nullptr) std::fclose(_file); }

	bool exists() const { return (_file != nullptr); }

	static const size_t IO_SIZE_MAX = 64 * 1024 * 1024;

	bool read(char * const ptr, const size_t size)
	{
		size_t s = size;
		char * cptr = ptr;
		while (s > IO_SIZE_MAX)
		{
			const size_t ret = std::fread(cptr, sizeof(char), IO_SIZE_MAX, _file);
			_crc32 = rc_crc32(_crc32, cptr, IO_SIZE_MAX);
			if (ret != IO_SIZE_MAX * sizeof(char))
			{
				error("failure of a read operation");
				return false;
			}
			cptr += IO_SIZE_MAX * sizeof(char);
			s -= IO_SIZE_MAX;
		}
		const size_t ret = std::fread(cptr, sizeof(char), s, _file);
		_crc32 = rc_crc32(_crc32, cptr, s);
		if (ret == s * sizeof(char)) return true;
		error("failure of a read operation");
		return false;
	}

	bool write(const char * const ptr, const size_t size)
	{
		size_t s = size;
		const char * cptr = ptr;
		while (s > IO_SIZE_MAX)
		{
			const size_t ret = std::fwrite(cptr, sizeof(char), IO_SIZE_MAX, _file);
			_crc32 = rc_crc32(_crc32, cptr, IO_SIZE_MAX);
			if (ret != IO_SIZE_MAX * sizeof(char))
			{
				error("failure of a write operation");
				return false;
			}
			cptr += IO_SIZE_MAX * sizeof(char);
			s -= IO_SIZE_MAX;
		}
		const size_t ret = std::fwrite(cptr, sizeof(char), s, _file);
		_crc32 = rc_crc32(_crc32, cptr, s);
		if (ret == s * sizeof(char)) return true;
		error("failure of a write operation");
		return false;
	}

	void write_crc32()
	{
		uint32_t crc32 = ~_crc32 ^ 0xa23777ac;
		write(reinterpret_cast<const char *>(&crc32), sizeof(crc32));
	}

	bool check_crc32()
	{
		uint32_t crc32 = 0, ocrc32 = ~_crc32 ^ 0xa23777ac;	// before the read operation
		read(reinterpret_cast<char *>(&crc32), sizeof(crc32));
		const bool success = (crc32 == ocrc32);
		if (!success) error("bad file (crc32)");
		return success;
	}
};
