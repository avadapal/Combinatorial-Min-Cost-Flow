#pragma once

#include <cassert>


template<class InternalQueueType>
class SSPVariantQueue {

public:
	using value_type = typename InternalQueueType::value_type;

	SSPVariantQueue(const InternalQueueType &s_in, const InternalQueueType &s_out) :
		s_in(s_in), s_out(s_out) {}

	SSPVariantQueue(InternalQueueType &&s_in, InternalQueueType &&s_out) :
		s_in(s_in), s_out(s_out) {}

	void push(const typename InternalQueueType::value_type &value, bool outgoing) {
		auto &queue = outgoing ? s_out : s_in;
		queue.push(value);
	}

	template<typename DeficitType>
	auto top(const DeficitType &deficit_s) -> typename InternalQueueType::value_type {
		assert(!empty());
		if (deficit_s < 0 || s_in.empty()) {
			assert(!s_out.empty());
			return s_out.top();
		} else {
			assert(!s_in.empty());
			return {-s_in.top().first, s_in.top().second};
		}
	}

	template<typename DeficitType>
	void pop(DeficitType deficit_s) {
		assert(!empty());
		auto &queue = deficit_s < 0 || s_in.empty() ? s_out : s_in;
		assert(!queue.empty());
		queue.pop();
	}
	
	auto empty() -> bool {
		return s_in.empty() && s_out.empty();
	}

private:
	InternalQueueType s_in;
	InternalQueueType s_out;

};


template<class InternalQueueType, typename Objective >
class SSPVariantAlternativeQueue {

public:
	using value_type = typename InternalQueueType::value_type;

 	SSPVariantAlternativeQueue( const Objective& obj ) : objective( obj ) {}
// 		s_in(s_in), s_out(s_out) {}
// 
// 	SSPVariantAlternativeQueue(InternalQueueType &&s_in, InternalQueueType &&s_out) :
// 		s_in(s_in), s_out(s_out) {}

	void push(const typename InternalQueueType::value_type &value, bool outgoing) {
		deltaS.push_back(value);
	}

	// TODO: this is not really a priority queue anymore since it just emits the elements from s_out first but doesn't compare s_in and s_out

	auto top() -> typename InternalQueueType::value_type {
		assert(!empty());
		auto first_less = [](const value_type& lhs, const value_type& rhs) {
				return lhs.first < rhs.first;
			};
		sort( deltaS.begin(), deltaS.end(), first_less );
 		cout << "|delta(S)| = " << deltaS.size() << endl;
		maximizer = deltaS.begin();
		auto max_delta = deltaS.front().first;
		auto max_value = objective( deltaS.front() );
//		cout << max_delta  << " " << max_value << endl;
		for ( auto it = maximizer+1, end = deltaS.end(); it != end; ++it ) {
			const auto &current_arc = *it;
			auto Delta = current_arc.first;
			auto next_value = objective( current_arc );
//			cout << Delta << " " << next_value << endl;
			if( next_value > max_value ) {
				max_value = next_value;
				max_delta = Delta;
				maximizer = it;
			}
		}
		
		cout << "maximum found: " << max_delta << " " << max_value << endl;

		return *maximizer;
	}

	void pop() {
		assert(!empty());
		deltaS.erase( maximizer );
	}

	auto empty() -> bool {
		return deltaS.empty();
	}

private:
	const Objective& objective;
	InternalQueueType deltaS;
  typename InternalQueueType::iterator maximizer;
};
