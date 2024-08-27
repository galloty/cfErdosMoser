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
#include <mutex>

#include <gmp.h>

// Memory allocation
class Heap
{
private:
	struct Stat
	{
		size_t size, max_size;
		size_t block_count, max_block_count, max_block_size;
		size_t alloc_count, realloc_count, free_count;

		Stat() : size(0), max_size(0), block_count(0), max_block_count(0), max_block_size(0), alloc_count(0), realloc_count(0), free_count(0) {}

		void reset()
		{
			max_size = size;
			max_block_count = block_count;
			max_block_size = 0;
			alloc_count = realloc_count = free_count = 0;
		}

		void alloc(const size_t alloc_size)
		{
			size += alloc_size;
			max_size = std::max(max_size, size);
			block_count += 1;
			max_block_count = std::max(max_block_count, block_count);
			max_block_size = std::max(max_block_size, alloc_size);
			alloc_count += 1;
		}

		void realloc(const size_t old_alloc_size, const size_t new_alloc_size)
		{
			size += new_alloc_size - old_alloc_size;
			max_size = std::max(max_size, size);
			max_block_size = std::max(max_block_size, new_alloc_size);
			realloc_count += 1;
		}

		void free(const size_t alloc_size)
		{
			size -= alloc_size;
			block_count -= 1;
			free_count += 1;
		}

	};

	static const size_t min_size = 256 * 1024 / sizeof(uint64_t);

	Stat _stat_g, _stat_ssg, _stat_gmp;
	std::queue<uint64_t *> _small_block_queue;
	std::queue<std::pair<uint64_t *, size_t>> _ssg_block_queue;
	std::queue<void *> _gmp_block_queue;
	std::mutex _mtx;

private:
	struct deleter { void operator()(const Heap * const p) { delete p; } };

public:
	Heap() { mp_set_memory_functions(_alloc_gmp, _realloc_gmp, _free_gmp); }
	virtual ~Heap() { mp_set_memory_functions(nullptr, nullptr, nullptr); }

	static Heap & get_instance()
	{
		static std::unique_ptr<Heap, deleter> p_instance(new Heap());
		return *p_instance;
	}

private:
	void _error(const std::string & who, const size_t size) const
	{
		std::cout << std::endl << get_memory_info1() << ", " << get_memory_info2() << ", " << get_memory_info3() << std::endl;
		std::ostringstream ss; ss << who << " failed, " << size << " bytes";
		throw std::runtime_error(ss.str());
	}

	static void * _aligned_alloc(const size_t size, const size_t alignment, const size_t offset = 0)
	{
		void * const alloc_ptr = std::malloc(size + alignment + offset + sizeof(size_t));
		if (alloc_ptr == nullptr) return nullptr;
		const size_t addr = reinterpret_cast<size_t>(alloc_ptr) + alignment + sizeof(size_t);
		size_t * const ptr = reinterpret_cast<size_t *>(addr - addr % alignment + offset);
		ptr[-1] = reinterpret_cast<size_t>(alloc_ptr);
		return (void *)(ptr);
	}

	static void _aligned_free(void * const ptr)
	{
		void * const alloc_ptr = reinterpret_cast<void *>(static_cast<size_t *>(ptr)[-1]);
		std::free(alloc_ptr);
	}

	uint64_t * _alloc(const size_t size)
	{
		const size_t alloc_size = size * sizeof(uint64_t);
		_stat_g.alloc(alloc_size);

		uint64_t * const ptr = static_cast<uint64_t *>(std::malloc(alloc_size));
		if (ptr == nullptr) _error("_alloc", alloc_size);
		return ptr;
	}

	uint64_t * _realloc(uint64_t * const ptr, const size_t old_size, const size_t new_size)
	{
		const size_t old_alloc_size = old_size * sizeof(uint64_t), new_alloc_size = new_size * sizeof(uint64_t);
		_stat_g.realloc(old_alloc_size, new_alloc_size);

		uint64_t * const new_ptr = static_cast<uint64_t *>(std::realloc(static_cast<void *>(ptr), new_alloc_size));
		if (new_ptr == nullptr) _error("_realloc", new_alloc_size);
		return new_ptr;
	}

	void _free(uint64_t * const ptr, const size_t size)
	{
		const size_t alloc_size = size * sizeof(uint64_t);
		_stat_g.free(alloc_size);

		std::free(static_cast<void *>(ptr));
	}

	uint64_t * _alloc_ssg(const size_t size)
	{
		const size_t alloc_size = size * sizeof(uint64_t);
		_stat_ssg.alloc(alloc_size);

		uint64_t * const ptr = static_cast<uint64_t *>(_aligned_alloc(alloc_size, 4096));	// 4kB TLB pages
		if (ptr == nullptr) _error("_alloc_ssg", alloc_size);
		return ptr;
	}

	void _free_ssg(uint64_t * const ptr, const size_t size)
	{
		const size_t alloc_size = size * sizeof(uint64_t);
		_stat_ssg.free(alloc_size);

		_aligned_free(static_cast<void *>(ptr));
	}

	static void * _alloc_gmp(size_t alloc_size)
	{
		auto & me = get_instance();
		std::lock_guard<std::mutex> lock(me._mtx);

		if (!me._gmp_block_queue.empty())
		{
			void * const ptr = me._gmp_block_queue.front(); me._gmp_block_queue.pop();
			size_t * const sptr = static_cast<size_t *>(ptr);
			const size_t size = sptr[0];
			if (size >= alloc_size) return static_cast<void *>(&sptr[1]);
			me._stat_gmp.free(size);
			std::free(ptr);
		}

		me._stat_gmp.alloc(alloc_size);
		void * const ptr = std::malloc(alloc_size + sizeof(size_t));
		if (ptr == nullptr) me._error("_alloc_gmp", alloc_size);
		size_t * const sptr = static_cast<size_t *>(ptr);
		sptr[0] = alloc_size;
		return static_cast<void *>(&sptr[1]);
	}

	static void * _realloc_gmp(void * old_ptr, size_t, size_t new_alloc_size)
	{
		auto & me = get_instance();
		std::lock_guard<std::mutex> lock(me._mtx);

		size_t * const old_sptr = static_cast<size_t *>(old_ptr);
		const size_t old_alloc_size = old_sptr[-1];

		if (old_alloc_size >= new_alloc_size) return old_ptr;

		me._stat_gmp.realloc(old_alloc_size, new_alloc_size);
		void * const new_ptr = std::realloc(static_cast<void *>(&old_sptr[-1]), new_alloc_size + sizeof(size_t));
		if (new_ptr == nullptr) me._error("_realloc_gmp", new_alloc_size);
		size_t * const new_sptr = static_cast<size_t *>(new_ptr);
		new_sptr[0] = new_alloc_size;
		return static_cast<void *>(&new_sptr[1]);
	}

	static void _free_gmp(void * ptr, size_t)
	{
		auto & me = get_instance();
		std::lock_guard<std::mutex> lock(me._mtx);

		size_t * const sptr = static_cast<size_t *>(ptr);
		me._gmp_block_queue.push(static_cast<void *>(&sptr[-1]));
	}

	static void get_unit(const size_t size, size_t & divisor, std::string & unit)
	{
		if (size < (size_t(100) << 10)) { divisor = 1; unit = "B"; }
		else if (size < (size_t(100) << 20)) { divisor = size_t(1) << 10; unit = "kB"; }
		else if (size < (size_t(100) << 30)) { divisor = size_t(1) << 20; unit = "MB"; }
		else { divisor = size_t(1) << 30; unit = "GB"; }
	}

	static size_t get_min_size(const size_t size) { return (size / Heap::min_size + 1) * Heap::min_size; }

public:
	static std::string get_size_str(const size_t size)
	{
		size_t size_divisor; std::string size_unit; get_unit(size, size_divisor, size_unit);
		std::ostringstream ss; ss << size / size_divisor << " " << size_unit;
		return ss.str();
	}

	std::string get_memory_size() const
	{
		std::ostringstream ss; ss << _stat_g.size << " + " << _stat_ssg.size << " + " << _stat_gmp.size << " B";
		return ss.str();
	}

	std::string get_memory_info1() const
	{
		const size_t size_g = _stat_g.max_size, size_ssg = _stat_ssg.max_size, size_gmp = _stat_gmp.max_size;
		const size_t block_size_g =  _stat_g.max_block_size, block_size_ssg = _stat_ssg.max_block_size, block_size_gmp = _stat_gmp.max_block_size;

		size_t size_divisor; std::string size_unit; get_unit(std::max(std::max(size_g, size_ssg), size_gmp), size_divisor, size_unit);
		size_t block_size_divisor; std::string block_size_unit; get_unit(std::max(std::max(block_size_g, block_size_ssg), block_size_gmp), block_size_divisor, block_size_unit);

		std::ostringstream ss;
		ss	<< "max size: " << size_g / size_divisor << " + " << size_ssg / size_divisor << " + " << size_gmp / size_divisor << " " << size_unit << ", "
			<< "max block size: " << block_size_g / block_size_divisor << " + " << block_size_ssg / block_size_divisor
			<< " + " << block_size_gmp / block_size_divisor << " " << block_size_unit;
		return ss.str();
	}

	std::string get_memory_info2() const
	{
		std::ostringstream ss;
		ss	<< "alloc: " << _stat_g.alloc_count << " + " << _stat_ssg.alloc_count << " + " << _stat_gmp.alloc_count << ", "
			<< "realloc: " << _stat_g.realloc_count << " + " << _stat_ssg.realloc_count << " + " << _stat_gmp.realloc_count << ", "
			<< "free: " << _stat_g.free_count << " + " << _stat_ssg.free_count << " + " << _stat_gmp.free_count;
		return ss.str();
	}

	std::string get_memory_info3() const
	{
		const size_t queue_size = _small_block_queue.size() * min_size;
		size_t size_divisor; std::string size_unit; get_unit(queue_size, size_divisor, size_unit);

		std::ostringstream ss;
		ss	<< "max block count: " << _stat_g.max_block_count << " + " << _stat_ssg.max_block_count << " + " << _stat_gmp.max_block_count << ", "
			<< "queue: " << _small_block_queue.size() << " (" << queue_size / size_divisor << " " << size_unit << ")";
		return ss.str();
	}

	std::string get_memory_usage() const
	{
		const size_t size_g = _stat_g.max_size, size_ssg = _stat_ssg.max_size, size_gmp = _stat_gmp.max_size;
		size_t size_divisor; std::string size_unit; get_unit(std::max(std::max(size_g, size_ssg), size_gmp), size_divisor, size_unit);
		std::ostringstream ss; ss << size_g / size_divisor << " + " << size_ssg / size_divisor << " + " << size_gmp / size_divisor << " " << size_unit << ")";
		return ss.str();
	}

	size_t get_max_mem_size() const { return _stat_g.max_size + _stat_ssg.max_size + _stat_gmp.max_size; }

	void reset()
	{
		while (!_small_block_queue.empty())
		{
			uint64_t * const ptr = _small_block_queue.front(); _small_block_queue.pop();
			_free(ptr, min_size);
		}

		while (!_ssg_block_queue.empty())
		{
			const auto blk = _ssg_block_queue.front(); _ssg_block_queue.pop();
			_free_ssg(blk.first, blk.second);
		}

		while (!_gmp_block_queue.empty())
		{
			void * const ptr = _gmp_block_queue.front(); _gmp_block_queue.pop();
			size_t * const sptr = static_cast<size_t *>(ptr);
			const size_t size = sptr[0];
			_stat_gmp.free(size);
			std::free(ptr);
		}

		_stat_g.reset(); _stat_ssg.reset(); _stat_gmp.reset();
	}

	uint64_t * alloc(size_t & size)
	{
		const size_t alloc_size = get_min_size(size);
		size = alloc_size;

		if ((alloc_size == min_size) && !_small_block_queue.empty())
		{
			uint64_t * const ptr = _small_block_queue.front(); _small_block_queue.pop();
			return ptr;
		}

		return _alloc(alloc_size);
	}

	uint64_t * realloc(uint64_t * const ptr, const size_t old_size, size_t & new_size, const bool copy)
	{
		if (old_size == new_size) return ptr;

		const size_t alloc_size = get_min_size(new_size);
		new_size = alloc_size;

		if ((alloc_size == min_size) && !_small_block_queue.empty())
		{
			uint64_t * const new_ptr = _small_block_queue.front(); _small_block_queue.pop();
			if (copy) for (size_t i = 0, n = std::min(alloc_size, old_size); i < n; ++i) new_ptr[i] = ptr[i];
			_free(ptr, old_size);	// because new_size = alloc_size = min_size != old_size
			return new_ptr;
		}

		if (old_size == min_size)
		{
			uint64_t * const new_ptr = _alloc(alloc_size);	// because new_size = alloc_size != old_size = min_size
			if (copy) for (size_t i = 0, n = std::min(alloc_size, old_size); i < n; ++i) new_ptr[i] = ptr[i];
			_small_block_queue.push(ptr);
			return new_ptr;
		}

		if (copy) return _realloc(ptr, old_size, alloc_size);
		_free(ptr, old_size);
		return _alloc(alloc_size);
	}

	void free(uint64_t * const ptr, const size_t size)
	{
		if (size == min_size) _small_block_queue.push(ptr);
		else _free(ptr, size);
	}

	uint64_t * alloc_ssg(size_t & size)
	{
		if (!_ssg_block_queue.empty())
		{
			const auto blk = _ssg_block_queue.front(); _ssg_block_queue.pop();
			if (blk.second >= size) { size = blk.second; return blk.first; }
			_free_ssg(blk.first, blk.second);
		}

		return _alloc_ssg(size);
	}

	void free_ssg(uint64_t * const ptr, const size_t size) { _ssg_block_queue.push(std::make_pair(ptr, size)); }
};
