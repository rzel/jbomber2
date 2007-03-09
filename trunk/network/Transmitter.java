package network;

import java.io.*;
import java.net.*;

import config.*;

/**
 * A <code>Transmitter</code> object sends data from a <code>Server</code> to a
 * connected client, it also recieves data from a connected client, and
 * forwards that data to the <code>Server</code>.
 */
public class Transmitter extends Thread
{
    // The socket that is connected to a client.
    private Socket client;

    // The server that created this Transmitter
    private Server server;

    // The InputStream to recieve data from the connected client.
    private ObjectInputStream obj_in;

    // The OutputStream to send data to the connected client.
    private ObjectOutputStream obj_out;

    // Stream to write log messages to
    private PrintStream log;

    // The id of the player to which this Transmitter is connected
    private int id;

    // Indicates whether this transmitter is recieving data
    private boolean running = true;

    /**
     * Creates a new <code>Transmitter</code> object and throws an
     * IOException when creating input- and outputstreams fails.
     *
     * @param client The <code>Socket</code> that is connected to a client.
     * @param server The <code>Server</code> that
     *               created this <code>Transmitter</code>.
     * @param log    <code>PrintStream</code> to write log messages to.
     * @param id     The unique id of the player that is connected to this
     *               <code>Transmitter</code>.
     * @throws IOException if an I/O error occurs when creating
     *                     input- and outputstreams.
     */
    public Transmitter(Socket client, Server server,
                       PrintStream log, int id) throws IOException
    {
        this.server = server;
        this.log    = log;
        this.client = client;
        this.id     = id;

        // Initialize connection
        OutputStream out = client.getOutputStream();
        obj_out = new ObjectOutputStream(out);

        InputStream in = client.getInputStream();
        obj_in = new ObjectInputStream(in);

        // Start recieving
        start();
    }

    /**
     * Recieves data from connected client.
     */
    public void run()
    {
        while (running)
        {
            try
            {
                Command command = (Command)obj_in.readObject();
                server.recieve(command);
                writeLog("Command recieved (" + command + ")");
                sleep(20);
            }
            catch (IOException ioe)
            {
                // Thrown when data could not be recieved.
                // Close this connection.
                close();
            }
            catch (ClassNotFoundException cnfe)
            {
                // Thrown when recieving data failed.
                // Close this connection.
                close();
            }
            catch (InterruptedException ie) {}
        }
    }

    /**
     * Sends data to a connected client. If sending fails,
     * this <code>Transmitter</code> is closed.
     */
    public void send(Command command)
    {
        try
        {
            obj_out.writeObject(command);
        }
        catch (IOException ioe)
        {
            writeLog("Error while sending, message" + ioe.getMessage());
            close();
        }
    }

    /**
     * Closes this connection and removes transmitter from server.
     * Also makes sure all other clients know this client left.
     */
    public void close()
    {
        running = false;

        // Make sure all other clients now this one left
        server.recieve(new PlayerCommand(id, Config.ACT_LEAVE));

        // Remove this client from the server
        server.removeClient(this);

        try
        {
            client.close();
            writeLog("Connection is closed");
        }
        catch (IOException ioe)
        {
             // Thrown if client was blocking on closing
            writeLog("Error while closing connection, message " +
                     ioe.getMessage());
        }
    }

    /**
     * Returns the closed state of the <code>Transmitter</code>.
     *
     * @Return <code>true</code> if the <code>Transmitter</code>
     *         has been closed.
     */
    public boolean isClosed()
    {
        return client.isClosed();
    }

    /**
     * Writes log messages to log if it exsits.
     *
     * @param message The message to write to log.
     */
    private void writeLog(String message)
    {
        if (log != null)
        {
            log.println(message);
        }
    }
}