#include "test.h"

namespace exst
{
    class IncidenceTests : public AbstractExstTest
    {
    CPPUNIT_TEST_SUITE(IncidenceTests);
            CPPUNIT_TEST(testSmallGraph);
            CPPUNIT_TEST(testIncidenceGraph);
        CPPUNIT_TEST_SUITE_END();
    private:
    public:
        void setUp()
        {
            AbstractExstTest::setUp();
        }

        void tearDown()
        {
            AbstractExstTest::tearDown();
        }

        void testSmallGraph()
        {
            std::list<lit_type> body;
            std::list<lit_type> head;
            body.push_back(*new lit_type(2, NEGATIVE));
            body.push_back(*new lit_type(3, NEGATIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count", (size_t) 2,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count", (size_t) 3,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            head.clear();
            body.clear();
            body.push_back(*new lit_type(4, NEGATIVE));
            body.push_back(*new lit_type(5, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count", (size_t) 4,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count", (size_t) 6,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            head.clear();
            body.clear();
            body.push_back(*new lit_type(5, POSITIVE));
            body.push_back(*new lit_type(6, POSITIVE));
            head.push_back(*new lit_type(7, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count", (size_t) 7,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count", (size_t) 9,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            head.clear();
            body.clear();
            head.push_back(*new lit_type(8, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count", (size_t) 8,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count", (size_t) 11,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            head.clear();
            body.clear();
            head.push_back(*new lit_type(9, NEGATIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count", (size_t) 9,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count", (size_t) 13,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            incidenceGraphStats->resetAssignment();

            (*atomIds)[(uint32_t) 1] = (uint32_t) 1;
            (*atomIds)[(uint32_t) 2] = (uint32_t) 2;
            (*atomIds)[(uint32_t) 3] = (uint32_t) 3;
            (*atomIds)[(uint32_t) 4] = (uint32_t) 4;
            (*atomIds)[(uint32_t) 5] = (uint32_t) 5;
            (*atomIds)[(uint32_t) 6] = (uint32_t) 6;
            (*atomIds)[(uint32_t) 7] = (uint32_t) 7;
            (*atomIds)[(uint32_t) 8] = (uint32_t) 8;
            (*atomIds)[(uint32_t) 9] = (uint32_t) 9;

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count Reduct", (size_t) 9,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count Reduct", (size_t) 13,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            lit_type l;
            l.id = 5;
            l.s = NEGATIVE;
            incidenceGraphStats->reduceGraph(l);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count Reduct", (size_t) 5,
                                         (incidenceGraphStats->getHTDIncidenceGraphReduct()->edgeCount()));
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count Reduct", (size_t) 13,
                                         incidenceGraphStats->getHTDIncidenceGraphReduct()->vertexCount());
        }

        void testIncidenceGraph()
        {
            std::list<lit_type> body;
            std::list<lit_type> head;
            body.push_back(*new lit_type(2, NEGATIVE));
            body.push_back(*new lit_type(3, NEGATIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            body.push_back(*new lit_type(4, NEGATIVE));
            body.push_back(*new lit_type(5, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            body.push_back(*new lit_type(5, POSITIVE));
            body.push_back(*new lit_type(6, POSITIVE));
            head.push_back(*new lit_type(7, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            head.push_back(*new lit_type(8, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            head.push_back(*new lit_type(9, NEGATIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            body.push_back(*new lit_type(6, POSITIVE));
            body.push_back(*new lit_type(4, POSITIVE));
            head.push_back(*new lit_type(2, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            body.push_back(*new lit_type(3, NEGATIVE));
            body.push_back(*new lit_type(5, POSITIVE));
            head.push_back(*new lit_type(7, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            body.push_back(*new lit_type(11, NEGATIVE));
            body.push_back(*new lit_type(9, NEGATIVE));
            head.push_back(*new lit_type(5, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            body.push_back(*new lit_type(6, POSITIVE));
            body.push_back(*new lit_type(15, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            head.clear();
            body.clear();
            body.push_back(*new lit_type(13, NEGATIVE));
            body.push_back(*new lit_type(9, POSITIVE));
            incidenceGraphStats->addRuleIncidenceGraph(body, head);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count", (size_t) 22,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count", (size_t) 21,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            incidenceGraphStats->resetAssignment();

            (*atomIds)[(uint32_t) 1] = (uint32_t) 1;
            (*atomIds)[(uint32_t) 2] = (uint32_t) 2;
            (*atomIds)[(uint32_t) 3] = (uint32_t) 3;
            (*atomIds)[(uint32_t) 4] = (uint32_t) 4;
            (*atomIds)[(uint32_t) 5] = (uint32_t) 5;
            (*atomIds)[(uint32_t) 6] = (uint32_t) 6;
            (*atomIds)[(uint32_t) 7] = (uint32_t) 7;
            (*atomIds)[(uint32_t) 8] = (uint32_t) 8;
            (*atomIds)[(uint32_t) 9] = (uint32_t) 9;
            (*atomIds)[(uint32_t) 10] = (uint32_t) 10;
            (*atomIds)[(uint32_t) 11] = (uint32_t) 11;
            (*atomIds)[(uint32_t) 12] = (uint32_t) 12;
            (*atomIds)[(uint32_t) 13] = (uint32_t) 13;
            (*atomIds)[(uint32_t) 14] = (uint32_t) 14;
            (*atomIds)[(uint32_t) 15] = (uint32_t) 15;

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count Reduct", (size_t) 22,
                                         incidenceGraphStats->getHTDIncidenceGraph()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count Reduct", (size_t) 21,
                                         incidenceGraphStats->getHTDIncidenceGraph()->vertexCount());

            lit_type l;
            l.id = 5;
            l.s = NEGATIVE;
            incidenceGraphStats->reduceGraph(l);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count Reduct", (size_t) 16,
                                         incidenceGraphStats->getHTDIncidenceGraphReduct()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count Reduct", (size_t) 21,
                                         incidenceGraphStats->getHTDIncidenceGraphReduct()->vertexCount());

            l.id = 3;
            l.s = POSITIVE;
            incidenceGraphStats->reduceGraph(l);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("Edge Count Reduct", (size_t) 14,
                                         incidenceGraphStats->getHTDIncidenceGraphReduct()->edgeCount());
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Node Count Reduct", (size_t) 21,
                                         incidenceGraphStats->getHTDIncidenceGraphReduct()->vertexCount());
        }
    };

    CPPUNIT_TEST_SUITE_REGISTRATION(IncidenceTests);
};