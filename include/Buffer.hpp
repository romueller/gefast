/*
 * SCT-PJ
 *
 * Copyright (C) 2016 Robert Mueller
 *
 * TODO add licence text (e.g. GNU (Affero) General Public License)
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#ifndef SCT_PJ_BUFFER_HPP
#define SCT_PJ_BUFFER_HPP

#include <atomic>
#include <mutex>
#include <queue>


namespace SCT_PJ {

/*
 * FIFO buffer suitable for concurrent accesses.
 */
template<typename T>
class Buffer {

public:
    Buffer() {
        flag_ = false;
    }

    void setFlag(const bool f) {
        flag_ = f;
    }

    bool getFlag() const {
        return flag_;
    }

    size_t size() {

        std::lock_guard<std::mutex> lock(mtx_);
        return buffer_.size();

    }

    T pop() {

        std::lock_guard<std::mutex> lock(mtx_);
        T val = buffer_.front();
        buffer_.pop();
        return val;

    }

    void push(T e) {

        std::lock_guard<std::mutex> lock(mtx_);
        buffer_.push(e);

    }


private:
    std::mutex mtx_;
    std::atomic<bool> flag_; // true signals that no more push-operations are expected

    std::queue<T> buffer_;

};


/*
 * Collection of 'Buffer' instances.
 *
 * The buffers take turns when new elements are pushed.
 * Closing the 'RotatingBuffers' instance signals to all underlying
 * buffers that further push-operations are not expected.
 */
template<typename T>
class RotatingBuffers {

public:
    RotatingBuffers(const unsigned long num = 1) {

        buffers_ = std::vector<Buffer<T>>(num);
        cur_ = 0;

    }

    void push(T t) {

        buffers_[cur_].push(t);
        cur_ = (cur_ + 1) % buffers_.size();

    }

    T pop(const unsigned long i) {

        return buffers_[i].pop();

    }

    Buffer<T>& getBuffer(const unsigned long i) {

        return buffers_[i];

    }

    void close() {

        for (auto iter = buffers_.begin(); iter != buffers_.end(); iter++) {
            iter->setFlag(true);
        }

    }


private:
    std::vector<Buffer<T>> buffers_;
    unsigned long cur_ = 0;

};

}

#endif //SCT_PJ_BUFFER_HPP