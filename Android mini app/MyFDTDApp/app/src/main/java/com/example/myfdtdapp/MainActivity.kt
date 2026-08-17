package com.example.myfdtdapp

import androidx.appcompat.app.AppCompatActivity
import android.os.Bundle
import android.widget.Button
import android.widget.TextView

class MainActivity : AppCompatActivity() {

    external fun runFDTD(): String

    companion object {
        init {
            System.loadLibrary("fdtd")
        }
    }

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)

        val btn = findViewById<Button>(R.id.buttonRun)
        val txt = findViewById<TextView>(R.id.textOutput)

        btn.setOnClickListener {
            txt.text = runFDTD()
        }
    }
}

