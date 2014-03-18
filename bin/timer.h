#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

/** @class Timer
 * @short Compute time between two points.
 * 
 * Use system time function.
 */
class Timer
{
public:
  Timer() { restart(); }
  void restart() //! Restart the timer
  {
    gettimeofday(&start_time, 0);
  }
   
  double s_elapsed()  //! Elapsed time in second
  {
    timeval end_time;
    gettimeofday(&end_time, 0);
    return end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec)/1000000.0;
  }
  
  double m_elapsed() //! Elapsed time in milli-second
  {
    timeval end_time;
    gettimeofday(&end_time, 0);
    return 1000.0*(end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec)/1000.0;
  }

  double u_elapsed() //! Elapsed time in micro-second
  {
    timeval end_time;
    gettimeofday(&end_time, 0);
    return 1000000.0*(end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec);
  }
  
  double elapsed() //! Elapsed time in milli-second
  {
    return m_elapsed();
  }

private:
  timeval start_time;
};

#endif // TIMER_H
