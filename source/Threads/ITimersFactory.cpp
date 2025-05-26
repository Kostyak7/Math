#include <Utils/Threads/ITimersFactory.h>

fem::ITimersFactory::TimerHandler::TimerHandler(std::weak_ptr<const ITimerHandlerImpl> impl)
    : m_impl(std::move(impl))
{
}

void fem::ITimersFactory::TimerHandler::stop_timer() {
    if (auto shared = m_impl.lock()) {
        shared->stop_timer_from_handler();
    }
}
