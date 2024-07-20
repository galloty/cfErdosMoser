/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>
#include <cstdlib> 
#include <string>
#include <sstream>
#include <queue>

#include <gmp.h>	// TODO remove

// Memory allocation
class Heap
{
private:
	static const size_t min_size = 64 * 1024 / sizeof(uint64_t);
	static size_t _size, _size_gmp;
	static size_t _alloc_count, _realloc_count, _free_count, _block_count;
	static size_t _max_size, _max_block_size, _max_block_count, _max_size_gmp, _max_block_size_gmp;
	static std::queue<uint64_t *> _small_block_queue;

	static uint64_t * _alloc(const size_t size)
	{
		++_alloc_count;
		uint64_t * const ptr = static_cast<uint64_t *>(std::malloc(size * sizeof(uint64_t)));
		if (ptr == nullptr) throw std::runtime_error("malloc failed");
		return ptr;
	}

	static uint64_t * _realloc(uint64_t * const ptr, const size_t size)
	{
		++_realloc_count;
		uint64_t * const new_ptr = static_cast<uint64_t *>(realloc(static_cast<void *>(ptr), size * sizeof(uint64_t)));
		if (new_ptr == nullptr) throw std::runtime_error("realloc failed");
		return new_ptr;
	}

	static void _free(uint64_t * const ptr)
	{
		++_free_count;
		free(static_cast<void *>(ptr));
	}

public:
	static size_t get_min_size(const size_t size) { return (size / Heap::min_size + 1) * Heap::min_size; }

	static uint64_t * alloc_function(const size_t size)
	{
		_size += size;

 		++_block_count;
		_max_size = std::max(_max_size, _size);
		_max_block_size = std::max(_max_block_size, size);
		_max_block_count = std::max(_max_block_count, _block_count);

		if ((size == min_size) && !_small_block_queue.empty())
		{
			uint64_t * const ptr = _small_block_queue.front(); _small_block_queue.pop();
			return ptr;
		}

		return _alloc(size);
	}

	static uint64_t * realloc_function(uint64_t * const ptr, const size_t old_size, const size_t new_size)
	{
		if (old_size == new_size) return ptr;

		_size += new_size - old_size;
		_max_size = std::max(_max_size, _size);
		_max_block_size = std::max(_max_block_size, new_size);

		if ((new_size == min_size) && !_small_block_queue.empty())
		{
			uint64_t * const new_ptr = _small_block_queue.front(); _small_block_queue.pop();
			for (size_t i = 0; i < min_size; ++i) new_ptr[i] = ptr[i];
			_free(ptr);
			return new_ptr;
		}

		if (old_size == min_size)
		{
			uint64_t * const new_ptr = _alloc(new_size);
			for (size_t i = 0; i < min_size; ++i) new_ptr[i] = ptr[i];
			_small_block_queue.push(ptr);
			return new_ptr;
		}

		return _realloc(ptr, new_size);
	}

	static void free_function(uint64_t * const ptr, const size_t size)
	{
		_size -= size;
		--_block_count;

		if (size == min_size) _small_block_queue.push(ptr);
		else _free(ptr);
	}

	static void * allocate_function_gmp(size_t size)
	{
		_size_gmp += size;
		_max_size_gmp = std::max(_max_size_gmp, _size_gmp);
		_max_block_size_gmp = std::max(_max_block_size_gmp, size);
		return malloc(size);
	}

	static void * reallocate_function_gmp(void *ptr, size_t old_size, size_t new_size)
	{
		_size_gmp += new_size - old_size;
		_max_size_gmp = std::max(_max_size_gmp, _size_gmp);
		_max_block_size_gmp = std::max(_max_block_size_gmp, new_size);
		return realloc(ptr, new_size);
	}

	static void free_function_gmp(void *ptr, size_t size) { _size_gmp -= size; if (size > 0) free(ptr); }

	static void get_unit(const size_t size, size_t & divisor, std::string & unit)
	{
		if (size < (size_t(10) << 10)) { divisor = 1; unit = "B"; }
		else if (size < (size_t(10) << 20)) { divisor = size_t(1) << 10; unit = "kB"; }
		else if (size < (size_t(10) << 30)) { divisor = size_t(1) << 20; unit = "MB"; }
		else { divisor = size_t(1) << 30; unit = "GB"; }
	}

public:
	Heap() { mp_set_memory_functions(allocate_function_gmp, reallocate_function_gmp, free_function_gmp); }
	virtual ~Heap() { mp_set_memory_functions(nullptr, nullptr, nullptr); }

	std::string	get_memory_size() const
	{
		std::ostringstream ss; ss << _size << " + " << _size_gmp << " B";
		return ss.str();
	}

	std::string	get_memory_info() const
	{
		const size_t max_size = _max_size * sizeof(uint64_t), max_size_gmp = _max_size_gmp;
		const size_t max_block_size = _max_block_size * sizeof(uint64_t), max_block_size_gmp = _max_block_size_gmp;

		size_t size_divisor; std::string size_unit; get_unit(std::max(max_size, max_size_gmp), size_divisor, size_unit);
		size_t block_size_divisor; std::string block_size_unit; get_unit(std::max(max_block_size, max_block_size_gmp), block_size_divisor, block_size_unit);

		std::ostringstream ss;
		ss << "max size: " << max_size / size_divisor << " + " << max_size_gmp / size_divisor << " " << size_unit << ", "
			<< "max block size: " << max_block_size / block_size_divisor << " + " << max_block_size_gmp / block_size_divisor << " " << block_size_unit << ", "
			<< "alloc: " << _alloc_count << ", realloc: " << _realloc_count << ", free: " << _free_count << ", max block count: " << _max_block_count
			<< " (" << _small_block_queue.size() << ").";
		return ss.str();
	}

	std::string	get_memory_usage() const
	{
		const size_t max_size = _max_size * sizeof(uint64_t), max_size_gmp = _max_size_gmp;
		size_t size_divisor; std::string size_unit; get_unit(std::max(max_size, max_size_gmp), size_divisor, size_unit);
		std::ostringstream ss; ss << max_size / size_divisor << " + " << max_size_gmp / size_divisor << " " << size_unit;
		return ss.str();
	}

	void reset()
	{
		_alloc_count = 0; _realloc_count = 0; _free_count = 0; _block_count = 0;
		_max_size = 0; _max_size_gmp = 0; _max_block_size = 0; _max_block_size_gmp = 0;

		while (!_small_block_queue.empty())
		{
			uint64_t * const ptr = _small_block_queue.front();
			_small_block_queue.pop();
			free(static_cast<void *>(ptr));
		}
	}
};

size_t Heap::_size = 0, Heap::_size_gmp = 0;
size_t Heap::_alloc_count = 0, Heap::_realloc_count = 0, Heap::_free_count = 0, Heap::_block_count = 0;
size_t Heap::_max_size = 0, Heap::_max_block_size = 0, Heap::_max_block_count = 0, Heap::_max_size_gmp = 0, Heap::_max_block_size_gmp = 0;
std::queue<uint64_t *> Heap::_small_block_queue;
