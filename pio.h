/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <thread>

class pio
{
private:
	struct deleter { void operator()(const pio * const p) { delete p; } };

public:
	pio() {}
	virtual ~pio() {}

	static pio & getInstance()
	{
		static std::unique_ptr<pio, deleter> pInstance(new pio());
		return *pInstance;
	}

private:
	void _print(const std::string & str) const
	{
		std::cout << str;
        std::ofstream logfile("cflog.txt", std::ios::app);
		if (logfile.is_open())
        {
            logfile << str;
            logfile.close();
        }
	}

public:
	static void print(const std::string & str) { getInstance()._print(str); }
};
