#pragma once

#include <queue>
#include <vector>
#include <utility>


template<
	class K,
	class V
>
class hybrid_queue {

public:
	using value_type = std::pair<K, V>;

	hybrid_queue():
		bucket(),
		pqueue([](const value_type &lhs, const value_type &rhs) { return lhs.first > rhs.first; }) {}

	void reserve_bucket(size_t new_cap) const {
		bucket.reserve(new_cap);
	}

	void push(const value_type &value) {
		const auto &key = value.first;
		if (key > min_key) {
			pqueue.push(value);
			return;
		}
		if (key == min_key) {
			bucket.push_back(value);
			return;
		}
		move_bucket_to_queue();
		bucket.push_back(value);
		min_key = key;
	}

	auto top() -> value_type {
		assert(!empty());
		if (bucket.empty())
			return pqueue.top();
#ifndef NDEBUG
		if (!pqueue.empty())
			for (const auto &entry : bucket)
				assert(entry.first <= pqueue.top().first);
#endif
		return bucket.back();
	}

	void pop() {
		assert(!empty());
		if (bucket.empty()) {
			min_key = pqueue.top().first;
			pqueue.pop();
		} else {
			bucket.pop_back();
		}
	}

	auto empty() const -> bool {
		return bucket.empty() && pqueue.empty();
	}

	void clear() {
		bucket.clear();
		pqueue.clear();
	}

private:

	void move_bucket_to_queue() {
		for (const auto &entry : bucket)
			pqueue.push(entry);
		bucket.clear();
	}

	K min_key;
	std::vector<value_type> bucket;
	std::priority_queue<value_type, std::vector<value_type>, bool (*)(const value_type &, const value_type &rhs)> pqueue;
};
