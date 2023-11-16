// ---------------------------------------------------------------------------------------------- //
// Includes
// ---------------------------------------------------------------------------------------------- //

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <optional>
#include <vector>
#include <algorithm>
#include <random>

// ---------------------------------------------------------------------------------------------- //
// Declarations
// ---------------------------------------------------------------------------------------------- //

// Namespace with types required to run the discrete event simulator.
namespace DES {

// An exponential random variable that can be sampled.
struct ExponentialRandomVariable { 
public:

    // Initializer
    ExponentialRandomVariable(double expectation); 

    // Expectation of the random variable.
    double expectation;

    // Generates and returns a random sample.
    double sample() const; 
private:

    // used when sampling
    static std::default_random_engine _engine;

    // used when sampling
    static std::uniform_real_distribution<double> _uniformRandomVariable;

};

// A packet of data that can be received or forwarded by a queue.
struct Packet { 
public:

    // Initializer
    Packet(int size, double timestamp); 
    
    // Size of the packet in bits.
    int size;

    // Timestamp of packet arrival or departure in seconds.
    double timestamp; 
};

// A line that is transmitting packets with specified average bit count and rate.
struct Line { 
public:

    // Initializer
    Line(double expectedSize, double expectedTransmissionRate); 

    // The amount of time the user has been receiving packets from the line.
    double runtime = 0;

    // Returns the expected size of the packets the line transmits
    double getExpectedSize() const;

    // Returns the expected rate packets are transmitted
    double getExpectedTransmissionRate() const;
    
    // Returns the packet at the receiving end of the line and update time.
    Packet getNextPacket(); 
private:

    // Random variable used to generate packet sizes.
    ExponentialRandomVariable _nextSize; 

    // Random variable used to generate packet arrival times.
    ExponentialRandomVariable _nextTime; 
};

// A data structure that represents a queue that can accept/forward packets from a line and output 
// overall performance reports.
struct Queue { 
public:

    // A report of performance measures taken by the queue.
    struct PerformanceReport { 
    public:

        // Initializer
        PerformanceReport(
            std::optional<int> capacity,
            double serviceRate,
            double runtime,
            double averageArrivalRate,
            double utilization, 
            double averageSize, 
            double idleProportion, 
            double dropProportion); 

        // The capacity of the queue in packets
        const std::optional<int> capacity;

        // The rate the queue forwards data in bits / s
        const double serviceRate;
        
        // Total time the queue has been running.
        const double runtime;

        // The rate the queue receives data in bits / s.
        const double averageArrivalRate;

        // the ratio between the arrival rate and service rate.
        const double utilization;

        // Average number of packets in the queue.
        const double averageSize; 

        // Proportion of time the queue is empty.
        const double idleProportion; 

        // Proportion of received packets that were dropped.
        const double dropProportion; 
    };

    // Initializer
    Queue(std::optional<int> capacity, double serviceRate);

    // Returns the capacity of the queue.
    std::optional<int> getCapacity() const;

    // Returns the service rate of the queue.
    double getServiceRate() const;

    // Returns a performance report of how the queue handled a line for a specified duration
    PerformanceReport connect(Line line, double duration); 
private:

    // Types of events that can occur on a queue
    enum EventType { arrival, departure, observer }; 

    // An event a queue can store in its history
    struct Event { 
    public:

        // Initializer
        Event(EventType type, double timestamp); 
        
        // Type of event.
        EventType type; 

        // Timestamp associated with event.
        double timestamp; 
    };

    // Number of packets the queue can store (std::nullopt ==> infinite size)
    std::optional<int> _capacity; 

    // Rate the queue forwards data in bits per second.
    double _serviceRate; 

    // Events that have occurred on the queue.
    std::vector<Event> _eventHistory = std::vector<Event>(); 

    // Total amount of time the queue has been running.
    double _runtime = 0; 

    // Time the last packet departed.
    double _lastDepartureTime = 0;

    // Total number of packets the queue has received.
    int _arrivalCount = 0; 

    // Total number of packets the queue has dropped.
    int _dropCount = 0; 

    // Total number of bits the queue has received.
    double _bitCount = 0;

    // Resets the queue to its initial state
    void _reset();

    // Adds a packet arrival and departure to event history.
    void _push(const Packet& packet); 

    // Returns a performance report of the queue at the current time.
    PerformanceReport _getPerformanceReport() const; 

    // Connects line to an infinite sized queue.
    void _connectInfinite(Line& line, double duration); 

    // Connects line to a finite sized queue.
    void _connectFinite(Line& line, double duration); 
};

// Writes a report buffer to a .csv file.
bool write(const std::vector<Queue::PerformanceReport>& reports, const std::string& filename);

} // namespace DES

// ---------------------------------------------------------------------------------------------- //
// Implementations
// ---------------------------------------------------------------------------------------------- //

namespace DES {

std::default_random_engine ExponentialRandomVariable::_engine = std::default_random_engine();

std::uniform_real_distribution<double> ExponentialRandomVariable::_uniformRandomVariable = std::uniform_real_distribution<double>(0,1);

ExponentialRandomVariable::ExponentialRandomVariable(double expectation):
expectation(expectation) {  }

double ExponentialRandomVariable::sample() const {
    return -expectation * log(1.0 - _uniformRandomVariable(_engine));
}

Packet::Packet(int size, double timestamp):
size(size),
timestamp(timestamp) {  }

Line::Line(double expectedSize, double expectedTransmissionRate):
_nextSize(ExponentialRandomVariable(expectedSize)),
_nextTime(ExponentialRandomVariable(1 / expectedTransmissionRate)) {  }

double Line::getExpectedSize() const {
    return _nextSize.expectation;
}

double Line::getExpectedTransmissionRate() const {
    return _nextTime.expectation;
}

Packet Line::getNextPacket() {
    runtime += _nextTime.sample();
    return Packet(int(_nextSize.sample()), runtime);
}

Queue::Queue(std::optional<int> capacity, double serviceRate):
_capacity(capacity),
_serviceRate(serviceRate) {  }

std::optional<int> Queue::getCapacity() const {
    return _capacity;
}

double Queue::getServiceRate() const {
    return _serviceRate;
}

Queue::PerformanceReport Queue::connect(Line line, double duration) {

    // reset
    line.runtime = 0;
    _reset();

    // determine how to connect line depending on queue size
    if(_capacity.has_value()) {
        _connectFinite(line, duration);
    } else {
        _connectInfinite(line, duration);
    }

    // get performance report
    const auto report = _getPerformanceReport();
    _reset();
    
    // return performance of connection
    return report;
}

void Queue::_reset() {
    _eventHistory.clear();
    _runtime = 0;
    _lastDepartureTime = 0;
    _arrivalCount = 0;
    _dropCount = 0;
    _bitCount = 0;
}

Queue::PerformanceReport Queue::_getPerformanceReport() const {

    // copy the event history
    auto events = _eventHistory;

    // remove departures that occur after runtime
    int i = int(events.size()) - 1;
    while(events[i].timestamp > _runtime) {
        events.erase(events.begin() + i);

        // arrivals and departures are interleaved
        // skip arrivals (which are known to occur during runtime)
        i -= 2;
    }
    const int departuresRemoved = (int(_eventHistory.size()) - (i + 1)) / 2;

    // add observers at 5 times the event frequency
    const double stepSize = _runtime / double(5 * events.size());
    for(double time = 0; time < _runtime; time += stepSize) {
        events.push_back(Event(EventType::observer, time));
    }
    const int observerCount = events.size() - _eventHistory.size() + departuresRemoved;

    // sort the events 
    std::sort(events.begin(), events.end(), [](const Event& event0, const Event& event1) {
        return event0.timestamp < event1.timestamp;
    });

    // run simulation
    const int eventsLastIndex = int(events.size()) - 1;
    int queueSize = 0;
    double queueSizeSum = 0;
    double timeIdle = 0;
    
    for(int i = 0; i < eventsLastIndex; ++i) {
        switch(events[i].type) {
        case EventType::arrival:
            queueSize += 1;
            break;
        case EventType::departure:
            queueSize -= 1;
            break;
        case EventType::observer:
            if(queueSize == 0) {
                timeIdle += events[i + 1].timestamp - events[i].timestamp;
            } else {
                queueSizeSum += queueSize;
            }
            break;
        }
    }

    // add last idle time if necessary
    if(events[eventsLastIndex].type == EventType::observer && queueSize == 0) {
        timeIdle += _runtime - events[eventsLastIndex].timestamp;
    }

    // return performance report
    return PerformanceReport(
        /* capacity: */ _capacity,
        /* serviceRate: */ _serviceRate,
        /* runtime: */ _runtime,
        /* averageArrivalRate: */ _bitCount / _runtime,
        /* utilization: */ _bitCount / _runtime / _serviceRate, 
        /* averageSize: */ queueSizeSum / double(observerCount), 
        /* idleProportion: */ timeIdle / _runtime, 
        /* dropProportion: */ double(_dropCount) / double(_arrivalCount));
}

void Queue::_push(const Packet& packet) {

    // update packet arrival count and bit count
    _arrivalCount += 1;
    _bitCount += packet.size;

    // add arrival event
    _eventHistory.push_back(Event(EventType::arrival, packet.timestamp));

    // add departure event
    double serviceTime = double(packet.size) / _serviceRate;
    if(packet.timestamp < _lastDepartureTime) {

        // packet received before last departure
        _lastDepartureTime +=  serviceTime;
        _eventHistory.push_back(Event(EventType::departure, _lastDepartureTime));
    } else {

        // packet received after last departure
        _lastDepartureTime = packet.timestamp + serviceTime;
        _eventHistory.push_back(Event(EventType::departure, _lastDepartureTime));
    }
}

void Queue::_connectInfinite(Line& line, double duration) {

    // receive the packets
    auto packet = line.getNextPacket();
    while(packet.timestamp < duration) {
        _push(packet);
        packet = line.getNextPacket();
    }

    // update runtime
    _runtime += duration;
}

void Queue::_connectFinite(Line& line, double duration) {

    // receive first packet
    auto packet = line.getNextPacket();
    _push(packet);

    // receive remaining packets
    const int queueCapacity = _capacity.value();
    int queueFront = 1; // stores the index of the departure time at the front
    int queueSize = 1; // current size of queue

    packet = line.getNextPacket();
    while(packet.timestamp < duration) {

        // remove packets departed since last event
        // note the arrivals and departures are interleaved in memory and chronological
        while(_eventHistory[queueFront].timestamp < packet.timestamp) {
            queueFront += 2; // move to next departure
            queueSize -= 1; // remove packet from queue

            if(queueFront >= _eventHistory.size()) {
                break;
            }
        }

        // accept packet if possible otherwise update drop count
        if(queueSize < queueCapacity) {
            queueSize += 1;
            _push(packet);
        } else {
            _arrivalCount += 1;
            _bitCount += packet.size;
            _dropCount += 1;
        }

        // get next packet
        packet = line.getNextPacket();
    }

    // update runtime
     _runtime += duration;
}

Queue::Event::Event(EventType type, double timestamp):
type(type),
timestamp(timestamp) {  }

Queue::PerformanceReport::PerformanceReport(
    std::optional<int> _capacity,
    double _serviceRate,
    double _runtime,
    double _averageArrivalRate,
    double _utilization, 
    double _averageSize, 
    double _idleProportion, 
    double _dropProportion):
capacity(_capacity),
serviceRate(_serviceRate),
runtime(_runtime),
averageArrivalRate(_averageArrivalRate),
utilization(_utilization),
averageSize(_averageSize),
idleProportion(_idleProportion),
dropProportion(_dropProportion) {  }

bool write(const std::vector<Queue::PerformanceReport>& reports, const std::string& filename) {

    // open file
    std::ofstream file;
    file.open(filename);
    if(!file.is_open()) {
        return false;
    }
    
    // write reports in .csv format
    file << "capacity,";
    file << "service_rate,";
    file << "runtime,";
    file << "average_arrival_rate,";
    file << "utilization,";
    file << "average_size,";
    file << "idle_proportion,";
    file << "drop_proportion" << std::endl;
    for(int i = 0; i < reports.size(); ++i) {
        if(reports[i].capacity.has_value()) {
            file << reports[i].capacity.value() << ",";
        } else {
            file << "infinity" << ",";
        }
        file << reports[i].serviceRate << ",";
        file << reports[i].runtime << ",";
        file << reports[i].averageArrivalRate << ",";
        file << reports[i].utilization << ",";
        file << reports[i].averageSize << ",";
        file << reports[i].idleProportion << ",";
        file << reports[i].dropProportion << std::endl;
    }
    
    // close file
    file.close();
    return true;
}

} // namespace DES

// ---------------------------------------------------------------------------------------------- //
// Simulations
// ---------------------------------------------------------------------------------------------- //

void simulation1() {
    const int sampleCount = 1000;
    const double lambda = 75;
    const double expectation = 1 / lambda;
    const double expectedVariance = expectation * expectation;

    std::cout << std::endl << "Simulation 1" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "sample count: " << sampleCount << std::endl;
    std::cout << "lambda: " << lambda << std::endl;
    std::cout << "expected average: " << expectation << std::endl;
    std::cout << "expected variance: " << expectedVariance << std::endl;

    // initialize samples
    auto randomVariable = DES::ExponentialRandomVariable(expectation);
    auto samples = std::vector<double>();
    for(int i = 0; i < sampleCount; ++i) {
        samples.push_back(randomVariable.sample());
    }
    
    // calculate and print average
    double average = 0;
    for(int i = 0; i < sampleCount; ++i) {
        average += samples[i];
    }
    average /= sampleCount;
    std::cout << std::endl << "average: " << average << std::endl;
    
    // calculate and print variance
    double variance = 0;
    for(int i = 0; i < sampleCount; ++i) {
        variance += pow(average - samples[i], 2);
    }
    std::cout << "variance: " << variance / sampleCount << std::endl;
}

void simulation2() {
    const double expectedPacketBitCount = 2000;
    const double serviceRate = 1000000;
    
    std::cout << std::endl << "Simulation 2" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    // make queue
    auto queue = DES::Queue(std::nullopt, serviceRate); // nullopt ==> infinite queue

    // run simulations
    auto reports = std::vector<DES::Queue::PerformanceReport>();

    for(double runtime = 1000; runtime < 3500; runtime += 1000) {
        for(double util = 0.25; util < 1.0; util += 0.1) {

            // make line
            double expectedTransmissionRate = util * serviceRate / expectedPacketBitCount;
            auto line = DES::Line(expectedPacketBitCount, expectedTransmissionRate);

            // connect line to queue and store performance report
            reports.push_back(queue.connect(line, runtime));
        }
    }
    
    // write results
    bool writeSucceeded = DES::write(reports, "q3_data.csv");
    if(writeSucceeded) {
        std::cout << "write succeeded" << std::endl;
    } else {
        std::cout << "write failed" << std::endl;
    }
}

void simulation3() {
    const double expectedPacketBitCount = 2000;
    const double serviceRate = 1000000;
    const double util = 1.2;
    const double expectedTransmissionRate = util * serviceRate / expectedPacketBitCount;
    
    std::cout << std::endl << "Simulation 3" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    // run simulation
    for(double runtime = 1000; runtime < 5500; runtime += 1000) {

        // create queue and line then generate report
        auto queue = DES::Queue(std::nullopt, serviceRate); // nullopt ==> infinite queue
        auto line = DES::Line(expectedPacketBitCount, expectedTransmissionRate);
        auto report = queue.connect(line, runtime);

        // output results
        std::cout << std::endl << "simulation time: " << runtime << std::endl;
        std::cout << "average queue size: " << report.averageSize << " packets" << std::endl;
        std::cout << "idle probability: " << report.idleProportion << std::endl;
    }
}

void simulation4() {
    const double expectedPacketBitCount = 2000;
    const double serviceRate = 1000000;

    std::cout << std::endl << "Simulation 4" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    
    auto queue = DES::Queue(10, serviceRate);
    auto reports = std::vector<DES::Queue::PerformanceReport>();

    for(double runtime = 1000; runtime < 3500; runtime += 1000) {
        for(int i = 0; i < 4; ++i) {

            if(i == 0) {
                queue = DES::Queue(10, serviceRate);
            } else if(i == 1) {
                queue = DES::Queue(25, serviceRate);
            } else if (i == 2) {
                queue = DES::Queue(50, serviceRate);
            } else if(i == 3) {
                queue = DES::Queue(75, serviceRate);
            }

            // run simulations
            for(double util = 0.5; util < 1.6; util += 0.1) {

                // make line
                double expectedTransmissionRate = serviceRate * util / expectedPacketBitCount;
                auto line = DES::Line(expectedPacketBitCount, expectedTransmissionRate);

                // connect line to queue and store performance report
                reports.push_back(queue.connect(line, runtime));
            }
        }
    }

    // write results
    bool write_succeeded = DES::write(reports, "q6_data.csv");
    if(write_succeeded) {
        std::cout << "write succeeded" << std::endl;
    } else {
        std::cout << "write failed" << std::endl;
    }
}

// ---------------------------------------------------------------------------------------------- //
// Program Entry
// ---------------------------------------------------------------------------------------------- //

int main(int argc, char* argv[]) {
    std::cout << std::endl << "Simulations" << std::endl;
    std::cout << "===============================================" << std::endl;

    simulation1();
    simulation2();
    simulation3();
    simulation4();
    
    std::cout << std::endl;
    std::cout << "===============================================" << std::endl;

    return 0;
}

