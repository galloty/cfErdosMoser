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

#include <gmp.h>

// Memory allocation
class Heap
{
private:
	static const size_t min_size = 64 * 1024 / sizeof(uint64_t);
	size_t _size, _size_other, _size_gmp;
	size_t _alloc_count, _realloc_count, _free_count, _block_count;
	size_t _max_size, _max_size_other, _max_size_gmp, _max_block_size, _max_block_count;
	std::queue<uint64_t *> _small_block_queue;

private:
	struct deleter { void operator()(const Heap * const p) { delete p; } };

public:
	Heap() : _size(0), _size_other(0), _size_gmp(0), _alloc_count(0), _realloc_count(0), _free_count(0), _block_count(0),
		_max_size(0), _max_size_other(0), _max_size_gmp(0), _max_block_size(0), _max_block_count(0)
	{
		mp_set_memory_functions(_alloc_gmp, _realloc_gmp, _free_gmp);
	}
	virtual ~Heap() { mp_set_memory_functions(nullptr, nullptr, nullptr); }

	static Heap & get_instance()
	{
		static std::unique_ptr<Heap, deleter> p_instance(new Heap());
		return *p_instance;
	}

private:
	static void * _aligned_alloc(const size_t size, const size_t alignment, const size_t offset = 0)
	{
		void * const alloc_ptr = std::malloc(size + alignment + offset + sizeof(size_t));
		const size_t addr = size_t(alloc_ptr) + alignment + sizeof(size_t);
		size_t * const ptr = (size_t *)(addr - addr % alignment + offset);
		ptr[-1] = size_t(alloc_ptr);
		return (void *)(ptr);
	}

	static void _aligned_free(void * const ptr)
	{
		void * const alloc_ptr = (void *)((size_t *)(ptr))[-1];
		std::free(alloc_ptr);
	}

	uint64_t * _alloc(const size_t size)
	{
		++_alloc_count;
		uint64_t * const ptr = static_cast<uint64_t *>(std::malloc(size * sizeof(uint64_t)));
		if (ptr == nullptr) throw std::runtime_error("malloc failed");
		return ptr;
	}

	uint64_t * _realloc(uint64_t * const ptr, const size_t size)
	{
		++_realloc_count;
		uint64_t * const new_ptr = static_cast<uint64_t *>(std::realloc(static_cast<void *>(ptr), size * sizeof(uint64_t)));
		if (new_ptr == nullptr) throw std::runtime_error("realloc failed");
		return new_ptr;
	}

	void _free(uint64_t * const ptr) { ++_free_count; std::free(static_cast<void *>(ptr)); }

	static void * _alloc_gmp(size_t size)
	{
		auto & me = get_instance();
		me._size_gmp += size;
		me._max_size_gmp = std::max(me._max_size_gmp, me._size_gmp);
		return std::malloc(size);
	}

	static void * _realloc_gmp(void * ptr, size_t old_size, size_t new_size)
	{
		auto & me = get_instance();
		me._size_gmp += new_size - old_size;
		me._max_size_gmp = std::max(me._max_size_gmp, me._size_gmp);
		return std::realloc(ptr, new_size);
	}

	static void _free_gmp(void * ptr, size_t size) { auto & me = get_instance(); me._size_gmp -= size; if (size > 0) std::free(ptr); }

	static void get_unit(const size_t size, size_t & divisor, std::string & unit)
	{
		if (size < (size_t(100) << 10)) { divisor = 1; unit = "B"; }
		else if (size < (size_t(100) << 20)) { divisor = size_t(1) << 10; unit = "kB"; }
		else if (size < (size_t(100) << 30)) { divisor = size_t(1) << 20; unit = "MB"; }
		else { divisor = size_t(1) << 30; unit = "GB"; }
	}

public:
	static std::string get_size_str(const size_t size)
	{
		size_t size_divisor; std::string size_unit; get_unit(size, size_divisor, size_unit);
		std::ostringstream ss; ss << size / size_divisor << " " << size_unit;
		return ss.str();
	}

	static size_t get_min_size(const size_t size) { return (size / Heap::min_size + 1) * Heap::min_size; }

	std::string get_memory_size() const
	{
		std::ostringstream ss; ss << _size << " + " << _size_other << " + " << _size_gmp << " B";
		return ss.str();
	}

	std::string get_memory_info() const
	{
		const size_t max_size = _max_size * sizeof(uint64_t), max_block_size = _max_block_size * sizeof(uint64_t);

		size_t size_divisor; std::string size_unit; get_unit(std::max(std::max(max_size, _max_size_other), _max_size_gmp), size_divisor, size_unit);
		size_t block_size_divisor; std::string block_size_unit; get_unit(max_block_size, block_size_divisor, block_size_unit);

		std::ostringstream ss;
		ss << "max size: " << max_size / size_divisor << " + " << _max_size_other / size_divisor << " + " << _max_size_gmp / size_divisor << " " << size_unit << ", "
			<< "max block size: " << max_block_size / block_size_divisor << " " << block_size_unit << ", "
			<< "alloc: " << _alloc_count << ", realloc: " << _realloc_count << ", free: " << _free_count << ", max block count: " << _max_block_count
			<< " (" << _small_block_queue.size() << ")";
		return ss.str();
	}

	std::string get_memory_usage() const
	{
		const size_t max_size = _max_size * sizeof(uint64_t);
		size_t size_divisor; std::string size_unit; get_unit(std::max(std::max(max_size, _max_size_other), _max_size_gmp), size_divisor, size_unit);
		std::ostringstream ss; ss << max_size / size_divisor << " + " << _max_size_other / size_divisor << " + " << _max_size_gmp / size_divisor << " " << size_unit;
		return ss.str();
	}

	size_t get_max_mem_size() const { return _max_size * sizeof(uint64_t) + _max_size_other + _max_size_gmp; }

	void reset()
	{
		_alloc_count = _realloc_count = _free_count = 0;
		_max_size = _max_size_other = _max_size_gmp = 0;
		_max_block_size = 0;
		_max_block_count = _block_count;

		while (!_small_block_queue.empty())
		{
			uint64_t * const ptr = _small_block_queue.front();
			_small_block_queue.pop();
			std::free(static_cast<void *>(ptr));
		}
	}

	uint64_t * alloc(const size_t size)
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

	uint64_t * realloc(uint64_t * const ptr, const size_t old_size, const size_t new_size)
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

	void free(uint64_t * const ptr, const size_t size)
	{
		_size -= size;
		--_block_count;

		if (size == min_size) _small_block_queue.push(ptr);
		else _free(ptr);
	}

	void * aligned_alloc(const size_t size)
	{
		_size_other += size;
		_max_size_other = std::max(_max_size_other, _size_other);

		void * const ptr = _aligned_alloc(size, 4096);	// 4kB TLB pages
		if (ptr == nullptr) throw std::runtime_error("malloc failed");
		return ptr;
	}

	void aligned_free(void * const ptr, const size_t size) { _size_other -= size; _aligned_free(ptr); }
};
