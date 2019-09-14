#include <iostream>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <signal.h>

struct Idle {
  void operator() () {
    long j = 0;
    for (long i = 0; i < 100000000; ++i) {
      j += i;
      if (i % 1000000 == 0) {
        std::cerr << j << std::endl;
      }
    }
    std::cerr << "Work done. Bye!" << std::endl;
    exit (0);
  }
};

struct Kill {
  template< typename T>
  void operator() (T t) {
    std::cerr << "Kill" << std::endl;
    exit (1);
  }
};


class timed_job
{
private:
  boost::asio::io_service io_service_;
  boost::asio::deadline_timer timer_;
  Idle idle;
  Kill kill;
  boost::thread t;

public:
  timed_job( int timeout ) :
    timer_( io_service_, boost::posix_time::seconds( timeout ) )  // Deadline timer
  {
  }

  void start()
  {
    timer_.async_wait (kill);

    // Post your work
    io_service_.post
      (
       boost::bind
       (
        &timed_job::do_work, this
        )
       );

    std::cout << "Not run yet." << std::endl;
    io_service_.run();
    std::cout << "stopped." << std::endl;
  }

private:
  void stop()
  {
    std::cerr << "call stop..." << std::endl;
    exit (1);
    io_service_.stop();
  }

  void do_work ()
  {

    std::cerr << "constructing..." << std::endl;

    t = boost::thread (idle);
    std::cerr << "constructed..." << std::endl;
    t.detach();
    std::cerr << "detached..." << std::endl;

    // // Keep posting the work.
    // io_service_.post
    //   (
    //    boost::bind
    //    (
    //     &timed_job::do_work, this
    //     )
    //    );
  }
};

int main()
{
  timed_job job( 1 );
  job.start();
  std::cout << "Normal exit" << std::endl;
  return 0;
}
