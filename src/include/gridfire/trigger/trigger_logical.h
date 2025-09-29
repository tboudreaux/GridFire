#pragma once

#include "gridfire/trigger/trigger_abstract.h"
#include "gridfire/trigger/trigger_result.h"

#include <string>
#include <vector>
#include <memory>

namespace gridfire::trigger {
    template <typename TriggerContextStruct>
    class LogicalTrigger : public Trigger<TriggerContextStruct> {};

    template <typename TriggerContextStruct>
    class AndTrigger final : public LogicalTrigger<TriggerContextStruct> {
    public:
        AndTrigger(std::unique_ptr<Trigger<TriggerContextStruct>> A, std::unique_ptr<Trigger<TriggerContextStruct>> B);
        ~AndTrigger() override = default;

        bool check(const TriggerContextStruct& ctx) const override;
        void update(const TriggerContextStruct& ctx) override;
        void reset() override;
        std::string name() const override;
        TriggerResult why(const TriggerContextStruct& ctx) const override;
        std::string describe() const override;
        size_t numTriggers() const override;
        size_t numMisses() const override;
    private:
        std::unique_ptr<Trigger<TriggerContextStruct>> m_A;
        std::unique_ptr<Trigger<TriggerContextStruct>> m_B;

        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;
    };

    template <typename TriggerContextStruct>
    class OrTrigger final : public LogicalTrigger<TriggerContextStruct> {
    public:
        OrTrigger(std::unique_ptr<Trigger<TriggerContextStruct>> A, std::unique_ptr<Trigger<TriggerContextStruct>> B);
        ~OrTrigger() override = default;

        bool check(const TriggerContextStruct& ctx) const override;
        void update(const TriggerContextStruct& ctx) override;
        void reset() override;
        std::string name() const override;
        TriggerResult why(const TriggerContextStruct& ctx) const override;
        std::string describe() const override;
        size_t numTriggers() const override;
        size_t numMisses() const override;

    private:
        std::unique_ptr<Trigger<TriggerContextStruct>> m_A;
        std::unique_ptr<Trigger<TriggerContextStruct>> m_B;

        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;
    };

    template <typename TriggerContextStruct>
    class NotTrigger final : public LogicalTrigger<TriggerContextStruct> {
    public:
        explicit NotTrigger(std::unique_ptr<Trigger<TriggerContextStruct>> A);
        ~NotTrigger() override = default;

        bool check(const TriggerContextStruct& ctx) const override;
        void update(const TriggerContextStruct& ctx) override;
        void reset() override;

        std::string name() const override;
        TriggerResult why(const TriggerContextStruct& ctx) const override;
        std::string describe() const override;
        size_t numTriggers() const override;
        size_t numMisses() const override;

    private:
        std::unique_ptr<Trigger<TriggerContextStruct>> m_A;

        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;

    };

    template <typename TriggerContextStruct>
    class EveryNthTrigger final : public LogicalTrigger<TriggerContextStruct> {
    public:
        explicit EveryNthTrigger(std::unique_ptr<Trigger<TriggerContextStruct>> A, size_t N);
        ~EveryNthTrigger() override = default;

        bool check(const TriggerContextStruct& ctx) const override;
        void update(const TriggerContextStruct& ctx) override;
        void reset() override;

        std::string name() const override;
        TriggerResult why(const TriggerContextStruct& ctx) const override;
        std::string describe() const override;
        size_t numTriggers() const override;
        size_t numMisses() const override;

    private:
        std::unique_ptr<Trigger<TriggerContextStruct>> m_A;

        size_t m_N;
        mutable size_t m_counter = 0;
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;
    };

    ///////////////////////////////
    // Templated Implementations //
    ///////////////////////////////


    template<typename TriggerContextStruct>
    AndTrigger<TriggerContextStruct>::AndTrigger(
        std::unique_ptr<Trigger<TriggerContextStruct>> A,
        std::unique_ptr<Trigger<TriggerContextStruct>> B
    ) : m_A(std::move(A)), m_B(std::move(B)) {}

    template<typename TriggerContextStruct>
    bool AndTrigger<TriggerContextStruct>::check(const TriggerContextStruct &ctx) const {
        const bool valid = m_A->check(ctx) && m_B->check(ctx);
        if (valid) {
            m_hits++;
        } else {
            m_misses++;
        }
        return valid;
    }

    template <typename TriggerContextStruct>
    void AndTrigger<TriggerContextStruct>::update(const TriggerContextStruct &ctx) {
        m_A->update(ctx);
        m_B->update(ctx);
        m_updates++;
    }

    template <typename TriggerContextStruct>
    void AndTrigger<TriggerContextStruct>::reset() {
        m_A->reset();
        m_B->reset();
        m_resets++;

        m_hits = 0;
        m_misses = 0;
        m_updates = 0;
    }

    template<typename TriggerContextStruct>
    std::string AndTrigger<TriggerContextStruct>::name() const {
        return "AND Trigger";
    }

    template<typename TriggerContextStruct>
    TriggerResult AndTrigger<TriggerContextStruct>::why(const TriggerContextStruct &ctx) const {
        TriggerResult result;

        result.name = name();

        TriggerResult A_result = m_A->why(ctx);
        result.causes.push_back(A_result);

        if (!A_result.value) {
            // Short Circuit
            result.value = false;
            result.description = "Failed because A (" + A_result.name + ") is false.";
            return result;
        }

        TriggerResult B_result = m_B->why(ctx);
        result.causes.push_back(B_result);

        if (!B_result.value) {
            result.value = false;
            result.description = "Failed because B (" + B_result.name + ") is false.";
            return result;
        }

        result.value = true;
        result.description = "Succeeded because both A (" + A_result.name + ") and B (" + B_result.description + ") are true.";
        return result;
    }

    template <typename TriggerContextStruct>
    std::string AndTrigger<TriggerContextStruct>::describe() const {
        return "(" + m_A->describe() + ") AND (" + m_B->describe() + ")";
    }

    template <typename TriggerContextStruct>
    size_t AndTrigger<TriggerContextStruct>::numTriggers() const {
        return m_hits;
    }

    template <typename TriggerContextStruct>
    size_t AndTrigger<TriggerContextStruct>::numMisses() const {
        return m_misses;
    }

    template <typename TriggerContextStruct>
    OrTrigger<TriggerContextStruct>::OrTrigger(
        std::unique_ptr<Trigger<TriggerContextStruct>> A,
        std::unique_ptr<Trigger<TriggerContextStruct>> B
    ) : m_A(std::move(A)), m_B(std::move(B)) {}

    template <typename TriggerContextStruct>
    bool OrTrigger<TriggerContextStruct>::check(const TriggerContextStruct &ctx) const {
        const bool valid = m_A->check(ctx) || m_B->check(ctx);
        if (valid) {
            m_hits++;
        } else {
            m_misses++;
        }
        return valid;
    }

    template <typename TriggerContextStruct>
    void OrTrigger<TriggerContextStruct>::update(const TriggerContextStruct &ctx) {
        m_A->update(ctx);
        m_B->update(ctx);
        m_updates++;
    }

    template <typename TriggerContextStruct>
    void OrTrigger<TriggerContextStruct>::reset() {
        m_A->reset();
        m_B->reset();
        m_resets++;

        m_hits = 0;
        m_misses = 0;
        m_updates = 0;
    }

    template<typename TriggerContextStruct>
    std::string OrTrigger<TriggerContextStruct>::name() const {
        return "OR Trigger";
    }

    template<typename TriggerContextStruct>
    TriggerResult OrTrigger<TriggerContextStruct>::why(const TriggerContextStruct &ctx) const {
        TriggerResult result;
        result.name = name();

        TriggerResult A_result = m_A->why(ctx);
        result.causes.push_back(A_result);

        if (A_result.value) {
            // Short Circuit
            result.value = true;
            result.description = "Succeeded because A (" + A_result.name + ") is true.";
            return result;
        }

        TriggerResult B_result = m_B->why(ctx);
        result.causes.push_back(B_result);

        if (B_result.value) {
            result.value = true;
            result.description = "Succeeded because B (" + B_result.name + ") is true.";
            return result;
        }

        result.value = false;
        result.description = "Failed because both A (" + A_result.name + ") and B (" + B_result.name + ") are false.";
        return result;
    }

    template <typename TriggerContextStruct>
    std::string OrTrigger<TriggerContextStruct>::describe() const {
        return "(" + m_A->describe() + ") OR (" + m_B->describe() + ")";
    }

    template <typename TriggerContextStruct>
    size_t OrTrigger<TriggerContextStruct>::numTriggers() const {
        return m_hits;
    }

    template <typename TriggerContextStruct>
    size_t OrTrigger<TriggerContextStruct>::numMisses() const {
        return m_misses;
    }

    template <typename TriggerContextStruct>
    NotTrigger<TriggerContextStruct>::NotTrigger(
        std::unique_ptr<Trigger<TriggerContextStruct>> A
    ) : m_A(std::move(A)) {}

    template <typename TriggerContextStruct>
    bool NotTrigger<TriggerContextStruct>::check(const TriggerContextStruct &ctx) const {
        const bool valid = !m_A->check(ctx);
        if (valid) {
            m_hits++;
        } else {
            m_misses++;
        }
        return valid;
    }

    template <typename TriggerContextStruct>
    void NotTrigger<TriggerContextStruct>::update(const TriggerContextStruct &ctx) {
        m_A->update(ctx);
        m_updates++;
    }

    template <typename TriggerContextStruct>
    void NotTrigger<TriggerContextStruct>::reset() {
        m_A->reset();
        m_resets++;

        m_hits = 0;
        m_misses = 0;
        m_updates = 0;
    }

    template<typename TriggerContextStruct>
    std::string NotTrigger<TriggerContextStruct>::name() const {
        return "NOT Trigger";
    }

    template<typename TriggerContextStruct>
    TriggerResult NotTrigger<TriggerContextStruct>::why(const TriggerContextStruct &ctx) const {
        TriggerResult result;
        result.name = name();

        TriggerResult A_result = m_A->why(ctx);
        result.causes.push_back(A_result);
        if (A_result.value) {
            result.value = false;
            result.description = "Failed because A (" + A_result.name + ") is true.";
            return result;
        }
        result.value = true;
        result.description = "Succeeded because A (" + A_result.name + ") is false.";
        return result;
    }

    template <typename TriggerContextStruct>
    std::string NotTrigger<TriggerContextStruct>::describe() const {
        return "NOT (" + m_A->describe() + ")";
    }

    template <typename TriggerContextStruct>
    size_t NotTrigger<TriggerContextStruct>::numTriggers() const {
        return m_hits;
    }

    template <typename TriggerContextStruct>
    size_t NotTrigger<TriggerContextStruct>::numMisses() const {
        return m_misses;
    }

    template <typename TriggerContextStruct>
    EveryNthTrigger<TriggerContextStruct>::EveryNthTrigger(std::unique_ptr<Trigger<TriggerContextStruct>> A, const size_t N) : m_A(std::move(A)), m_N(N) {
        if (N == 0) {
            throw std::invalid_argument("N must be greater than 0.");
        }
    }

    template <typename TriggerContextStruct>
    bool EveryNthTrigger<TriggerContextStruct>::check(const TriggerContextStruct &ctx) const
    {
        if (m_A->check(ctx) && (m_counter % m_N == 0)) {
            m_hits++;
            return true;
        }
        m_misses++;
        return false;
    }

    template <typename TriggerContextStruct>
    void EveryNthTrigger<TriggerContextStruct>::update(const TriggerContextStruct &ctx) {
        if (m_A->check(ctx)) {
            m_counter++;
        }
        m_A->update(ctx);
        m_updates++;
    }

    template <typename TriggerContextStruct>
    void EveryNthTrigger<TriggerContextStruct>::reset() {
        m_A->reset();
        m_resets++;

        m_counter = 0;
        m_hits = 0;
        m_misses = 0;
        m_updates = 0;
    }

    template <typename TriggerContextStruct>
    std::string EveryNthTrigger<TriggerContextStruct>::name() const {
        return "Every Nth Trigger";
    }

    template <typename TriggerContextStruct>
    TriggerResult EveryNthTrigger<TriggerContextStruct>::why(const TriggerContextStruct &ctx) const {
        TriggerResult result;
        result.name = name();
        TriggerResult A_result = m_A->why(ctx);
        result.causes.push_back(A_result);
        if (!A_result.value) {
            result.value = false;
            result.description = "Failed because A (" + A_result.name + ") is false.";
            return result;
        }
        if (m_counter % m_N == 0) {
            result.value = true;
            result.description = "Succeeded because A (" + A_result.name + ") is true and the counter (" + std::to_string(m_counter) + ") is a multiple of N (" + std::to_string(m_N) + ").";
            return result;
        }
        result.value = false;
        result.description = "Failed because the counter (" + std::to_string(m_counter) + ") is not a multiple of N (" + std::to_string(m_N) + ").";
        return result;
    }

    template <typename TriggerContextStruct>
    std::string EveryNthTrigger<TriggerContextStruct>::describe() const {
        return "Every " + std::to_string(m_N) + "th (" + m_A->describe() + ")";
    }

    template <typename TriggerContextStruct>
    size_t EveryNthTrigger<TriggerContextStruct>::numTriggers() const {
        return m_hits;
    }

    template <typename TriggerContextStruct>
    size_t EveryNthTrigger<TriggerContextStruct>::numMisses() const {
        return m_misses;
    }



}
