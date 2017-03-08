/*
 * SCT-PJ
 *
 * Copyright (C) 2016 - 2017 Robert Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#ifndef SCT_PJ_BUFFER_HPP
#define SCT_PJ_BUFFER_HPP

#include <condition_variable>
#include <mutex>
#include <queue>


namespace SCT_PJ {

/*
 * FIFO buffer (also) suitable for concurrent accesses.
 */
template<typename T>
class Buffer {

public:
    Buffer() {
        closed_ = false;
    }

    inline void close() {
        closed_ = true;
    }

    inline bool isClosed() const {
        return closed_;
    }

    inline size_t size() {
        return buffer_.size();
    }

    inline T pop() {

        T val = buffer_.front();
        buffer_.pop();
        return val;

    }

    inline void push(T& t) {
        buffer_.push(t);
    }

    inline void push(std::vector<T>& v) {

        for (auto iter = v.begin(); iter != v.end(); iter++) {
            buffer_.push(*iter);
        }

    }

    inline void swapContents(Buffer<T>& b) {
        buffer_.swap(b.buffer_);
    }

    inline void syncClose() {

        std::unique_lock<std::mutex> lock(mtx_);
        close();
        lock.unlock();
        cv_.notify_one();

    }

    inline bool syncIsClosed() {

        std::unique_lock<std::mutex> lock(mtx_);
        return isClosed();

    }

    inline size_t syncSize() {

        std::unique_lock<std::mutex> lock(mtx_);
        return size();

    }

    inline T syncPop() {

        std::unique_lock<std::mutex> lock(mtx_);
//        cv_.wait(lock, [this](){return buffer_.size() > 0 || flag_;});
        return pop();

    }

    inline void syncPush(T& t) {

        std::unique_lock<std::mutex> lock(mtx_);
        push(t);
        lock.unlock();
        cv_.notify_one();

    }

    inline void syncPush(std::vector<T>& v) {

        std::unique_lock<std::mutex> lock(mtx_);
        push(v);
        lock.unlock();
        cv_.notify_one();

    }

    inline void syncSwapContents(Buffer<T>& b) {//assumption: no concurrent accesses on b (i.e. b is a (thread-)local buffer)

        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [this](){return buffer_.size() > 0 || closed_;});
        swapContents(b);

    }

private:
    std::mutex mtx_;
    std::condition_variable cv_;
    bool closed_; // true signals that no more push-operations are expected

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

    void push(T& t) {

        buffers_[cur_].syncPush(t);
        cur_ = (cur_ + 1) % buffers_.size();

    }

    void push(std::vector<T>& t) {

        buffers_[cur_].syncPush(t);
        cur_ = (cur_ + 1) % buffers_.size();

    }

    T pop(const unsigned long i) {
        return buffers_[i].syncPop();
    }

    Buffer<T>& getBuffer(const unsigned long i) {
        return buffers_[i];
    }

    void close() {

        for (auto iter = buffers_.begin(); iter != buffers_.end(); iter++) {
            iter->syncClose();
        }

    }


private:
    std::vector<Buffer<T>> buffers_;
    unsigned long cur_;

};

}

#endif //SCT_PJ_BUFFER_HPP