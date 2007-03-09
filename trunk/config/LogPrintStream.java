package config;

import java.io.PrintStream;
import java.io.OutputStream;

import javax.swing.JTextArea;

/**
 * A <code>LogPrintStream</code> object writes data
 * to a <code>JTextArea</code>.
 * Currently only the <code>println()</code> methods are overriden.
 */
public class LogPrintStream extends PrintStream
{
    // The JTextArea to print messages to
    private JTextArea area;

    /**
     * Creates a new <code>LogPrintStream</code> object.
     *
     * @param area The <code>JTextArea</code> to print to.
     */
    public LogPrintStream(JTextArea area)
    {
        super(System.out);
        this.area = area;
    }

    /**
     * Print a message to the JTextArea, followed by a new line.
     *
     * @param message The message to print.
     */
    public void println(String message)
    {
        area.append(message);
        println();
    }

    /**
     * Print a message new lint JTextArea.
     */
    public void println()
    {
        area.append("\n");
    }
}