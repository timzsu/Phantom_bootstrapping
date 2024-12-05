#pragma once

#include <chrono>
#include <iostream>

namespace nexus {

class Timer {
 private:
  std::chrono::_V2::system_clock::time_point start_;
  std::chrono::_V2::system_clock::time_point end_;

  template <typename T = std::chrono::milliseconds>
  inline void print_duration(const char *text) {
    std::cout << text << std::chrono::duration_cast<T>(end_ - start_).count() << std::endl;
  }

 public:
  Timer() {
    start_ = std::chrono::high_resolution_clock::now();
  }

  inline void start() {
    start_ = std::chrono::high_resolution_clock::now();
  }

  template <typename T = std::chrono::milliseconds>
  inline void stop(const char *text) {
    end_ = std::chrono::high_resolution_clock::now();
    print_duration<T>(text);
  }

  inline void stop() {
    end_ = std::chrono::high_resolution_clock::now();
  }

  template <typename T = std::chrono::milliseconds>
  inline long duration() {
    return std::chrono::duration_cast<T>(end_ - start_).count() / 1.0;
  }
};
}  // namespace nexus
