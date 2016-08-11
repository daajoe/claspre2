#include <exst/ExstTypes.h>
#include "test.h"

namespace exst
{
    class ExstTypesTests : public AbstractExstTest
    {
    CPPUNIT_TEST_SUITE(ExstTypesTests);
            CPPUNIT_TEST(testCopyGraph);
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

        void testCopyGraph()
        {
            MyGraph graph, graphCopy;
            std::unordered_map<uint32_t, EdgeType> node1;
            std::unordered_map<uint32_t, EdgeType> node2;
            node1[2] = exst::HEAD;
            node2[1] = exst::BODY;
            graph[1] = node1;
            graph[2] = node2;
            graphCopy = copyMyGraph(graph);
            CPPUNIT_ASSERT_EQUAL_MESSAGE("Assert size of copied graph is equal to original", graph.size(),
                                         graphCopy.size());
        }
    };

    CPPUNIT_TEST_SUITE_REGISTRATION(ExstTypesTests);
};