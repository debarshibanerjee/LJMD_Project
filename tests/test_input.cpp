#include "prototypes.h"
#include "gtest/gtest.h"

TEST(test_input, get_a_line) {
	char line[BLEN], buffer[BLEN];
	FILE* fstream;

	sprintf(buffer, "100\n");
	fstream = fmemopen(buffer, BLEN, "r");
	ASSERT_NE(fstream, nullptr);

	get_a_line(fstream, line);
	ASSERT_EQ(100, atoi(line));
	fclose(fstream);

	sprintf(buffer, "101.357\n");
	fstream = fmemopen(buffer, BLEN, "r");
	ASSERT_NE(fstream, nullptr);

	get_a_line(fstream, line);
	ASSERT_DOUBLE_EQ(101.357, atof(line));
	fclose(fstream);

	sprintf(buffer, "random string\n");
	fstream = fmemopen(buffer, BLEN, "r");
	ASSERT_NE(fstream, nullptr);

	get_a_line(fstream, line);
	ASSERT_STREQ("random string", line);
	fclose(fstream);
}

